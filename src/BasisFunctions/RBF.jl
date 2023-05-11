"""
    RBF{IT<:Integer , FT<:AbstractFloat} <: LinearBasisFunction

屋顶基函数 (Rooftop basis function, RBF) 基函数复合类型：
```
isbd        ::Bool              是否为边界元即半基函数，布尔类型
bfID        ::IT                基函数编号，整形
inGeo       ::MVector{2, IT}    基函数所在两个六面体编号（半基函数为1个，赋值成一样的两个），长度为2的向量数组
inGeoID     ::MVector{2, IT}    基函数在两个六面体中的局部编号（半基函数为1个，赋值成一样的两个），取值1:4，长度为2的向量数组
center      ::MVec3D{FT}        基函数中心，用于八叉树分组
```
"""
mutable struct RBF{IT<:Integer , FT<:AbstractFloat} <: LinearBasisFunction
    isbd        ::Bool
    bfID        ::IT
    inGeo       ::MVector{2, IT}
    inGeoID     ::MVector{2, IT}
    center      ::MVec3D{FT}
end

"""
    RBF{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}

RBF 的默认构造函数，所有元素置零。
"""
function RBF{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    isbd        =    true
    bfID        =    zero(IT)
    inGeo       =    zero(MVector{2, IT})
    inGeoID     =    zero(MVector{2, IT})
    center      =    zero(MVec3D{FT})
    return RBF{IT, FT}(isbd, bfID, inGeo, inGeoID, center)
end
RBF()   =    RBF{IntDtype, Precision.FT}()

@doc raw"""
    boundaryRBF(h::Hexahedron)

重载面的提取顺序以匹配屋顶基函数 (RBF) 在六面体中的面按 
```math
u=1, u=0, v = 1, v = 0, w = 1, w = 0 
```
的顺序排列。
"""
function boundaryRBF(h::Hexahedron)
    indices = [ (2,3,7,6),(1,4,8,5),(4,3,7,8),
                (1,2,6,5),(5,6,7,8),(1,2,3,4)]
    SimpleMesh(h.vertices, connect.(indices))
end


"""
    setQuad4Hexas!(hexameshData::HexahedraMesh{IT, FT}, hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}, ::Val{:RBF}) where {IT, FT, CT}

计算构成六面体的所有四边形，并写入六面体 `hexasInfo`，给屋顶基函数 (RBF) 基函数赋值。
"""
function setQuad4Hexas!(hexameshData::HexahedraMesh{IT, FT}, hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}, ::Val{:RBF}) where {IT, FT, CT}
    # 六面体数
    hexnum      =   hexameshData.hexnum
    # 组成每个六面体的点的 id 数组
    hexasIDs    =   hexameshData.hexahedras
    # 组成每个面的局部 id 
    #* 注意不可改变下面的顺序因为跟 RBF 基函数对应的面在六面体的中顺序有关 *#
    quad6LocalIndices = [[2,3,7,6],[1,4,8,5],[4,3,7,8],[1,2,6,5],[5,6,7,8],[1,2,3,4]]

    #=  创建数组，用于后续保存所有面的点信息(4*2, 前四列排序
        寻找单独面，后四列保持四边形的点顺序)，然后是面所在的六面体id(1)，
        面所在 uniquad(独立四边形) 编号(1)、正负等信息(2)，共 12 列      =#
    face4hexd2bf2   =   zeros(IT, (hexnum*6, 12))
    # 利用内存视图提取其中的所有面的点
    face4           =   reshape(view(face4hexd2bf2, :, 1:4), (:, 24))
    face4quadids    =   reshape(view(face4hexd2bf2, :, 5:8), (:, 24))
    # 利用内存视图提取六面体信息的tri2和基函数信息的unitri2
    hex2            =   @view face4hexd2bf2[:, 9:10]
    uniquad2        =   @view face4hexd2bf2[:, 11:12]

    ## 创建六面体索引( id+localid )并根据排序更新
    # 此维度为六面体id
    hex2[:, 1]  .=   repeat(1:hexnum,  outer = 6 )
    # 此维度为该面在所在六面体的局部id
    hex2[:, 2]  .=   repeat(1:6,  inner = hexnum )
    
    # 将每个六面体的六个个面的点存放在 face4 数组中
    @threads for facei in 1:6
        viewofnodeids = view(hexasIDs, quad6LocalIndices[facei], :)'
        face4[:, facei:6:end] .=    viewofnodeids
        face4quadids[:, facei:6:end] .= viewofnodeids
    end
    # 重塑形状将每个边点组合放入同一维度
    face4        =  reshape(face4, (:, 4))
    face4quadids =  reshape(face4quadids, (:, 4))

    ##  对前四排序使得同一面的表示形式相同
    #   对前四个构成面的点排序
    @views sort!(face4hexd2bf2[:, 1:4], dims = 2)
    # 将每个六面体的四个角标点作为整体排序，次序按v1,v2,v3
    face4hexd2bf2   .=  sortslices(face4hexd2bf2, dims = 1, by = x -> (x[1], x[2], x[3], x[4]), alg = ThreadsX.MergeSort)

    ## for 循环计算PWC基函数id
    @inbounds let uniquadplusid   =   1
        uniquad2[1, 1]  =   1; uniquad2[1, 2] = 1
        # 初始值即1的位置要预计算一下
        uniquadplusid   +=  1

        # 从2开始循环，跟前一边不同为 unitri+，跟前一边不同为 unitri-
        for ii in 2:(hexnum*6)
            if  (face4[ii, 1]  != face4[ii-1, 1]) ||
                (face4[ii, 2]  != face4[ii-1, 2]) ||
                (face4[ii, 3]  != face4[ii-1, 3]) ||
                (face4[ii, 4]  != face4[ii-1, 4])
                uniquad2[ii, 1]      =   uniquadplusid
                uniquad2[ii, 2]      =   1
                uniquadplusid       +=   1
            elseif  (face4[ii, 1]  ==   face4[ii-1, 1]) &&
                (face4[ii, 2]  ==   face4[ii-1, 2]) && 
                (face4[ii, 3]  ==   face4[ii-1, 3]) && 
                (face4[ii, 4]  ==   face4[ii-1, 4]) 
                uniquad2[ii, 1]      =   uniquadplusid  - 1
                uniquad2[ii, 2]      =   -1
            end #if
        end #for
    end #begin

    # uniquadnum 函数数量( 同时也是独立四边形的数量 )
    uniquadnum   =   uniquad2[end, 1]
    # 单独四边形所在的 face4 的行 id
    uniquad2nodeRow    =   zeros(IT, uniquadnum)
    iunitri     =   1
    for inode in 1:(hexnum*6)
        (uniquad2[inode, 2] == 1) && begin
            uniquad2nodeRow[iunitri]     =   inode
            iunitri += 1
        end
    end
    # 预分配所有的构成六面体的四边形
    uniqueQuads  =   [Quads4Hexa{IT, FT}()  for _ in 1:uniquadnum ]

    # 设置 uniqueQuads 的三个点、边长、边法向量等信息
    @views for iQuad in 1:uniquadnum
        # 该四边形
        quad    =   uniqueQuads[iQuad]
        vertices    = quad.vertices

        # 该四边形对应的四个角点id
        quadvertIDs  =   view(face4quadids, uniquad2nodeRow[iQuad], :)
        # 设置角点
        for iv in 1:4
            vertices[:, iv]     .=   view(hexameshData.node, :, quadvertIDs[iv])
        end

        # 边指向向量
        edgev̂          =   quad.edgev̂
        @views for edgei in 1:4
            edgev̂[:, edgei]    =   vertices[:, edgei%4+1]  - vertices[:, edgei]
        end

        # 边长
        for edgei in 1:4
            @views quad.edgel[edgei]   =   norm(edgev̂[:, edgei])
        end

        @views facen̂    =  -cross(edgev̂[:, 1], edgev̂[:, 2])
        facen̂    ./=   norm(facen̂)
        # 计算边的指向单位向量、外法单位向量
        edgen̂                  =   quad.edgen̂
        for edgei in 1:4
            # 边的指向向量单位化
            @views edgev̂[:, edgei]  ./=   quad.edgel[edgei]
            # 最后计算边外法向量
            @views edgen̂[:, edgei]    =   cross(facen̂, edgev̂[:, edgei])
        end
    end

    ## 开始构造 RBF 基函数
    # rbf 函数数量
    rbfnum      =   uniquad2[end, 1]
    # 保存 rbf 基函数信息的数组
    rbfsInfo     =   [RBF{IT, FT}() for _ in 1:rbfnum]
    # 写入公共面数是否为边界元信息、基函数相关信息
    for faceii in 1:6hexnum
        # 面id
        faceID      =   uniquad2[faceii, 1]
        # 正负？
        ispOrm      =   uniquad2[faceii, 2]
        # 对应基函数
        rbfInfo     =   rbfsInfo[faceID]
        # 对应面（四边形）
        uniqueQaud  =   uniqueQuads[faceID]
        # 基函数序号
        rbfInfo.bfID    =   faceID
        # 基函数中心
        rbfInfo.center[:]   =   mean(uniqueQaud.vertices, dims = 2)
        # 其它信息
        (ispOrm == 1) && begin
            #  + 基函数所在六面体、是该六面体第几个基函数
            rbfInfo.inGeo[1]        =   hex2[faceii,1]
            rbfInfo.inGeoID[1]      =   hex2[faceii,2]
        end
        (ispOrm == -1) && begin
            # 有 - 基函数的必为非边界元
            uniqueQaud.isbd     =   false
            rbfInfo.isbd        =   false
            #  - 基函数所在六面体、是该六面体第几个基函数
            rbfInfo.inGeo[2]        =   hex2[faceii,1]
            rbfInfo.inGeoID[2]      =   hex2[faceii,2]
        end
    end
    ## 设置六面体所包含的基函数信息
    # 按六面体编号排列，效果是相邻的四个四边形同属一个六面体并桉顺序排列
    face4hexd2bf2   .=  sortslices(face4hexd2bf2, dims = 1, by = x -> (x[9], x[10]), alg = ThreadsX.MergeSort)

    rbfsID  =   reshape(uniquad2[:, 1], (6, :))
    # 找出每个六面体对应的六个面，并在 “-” 面将四边形面积变负
    for iHexa in eachindex(hexasInfo)
        # 该六面体
        hexa    =   hexasInfo[iHexa]
        # 四个四边形 id 
        hexa.faces      .=   uniqueQuads[uniquad2[(6(iHexa -1)+1):(6iHexa), 1]]
        # 将在负半部分的面积取反
        hexa.facesArea .*=   @view uniquad2[(6(iHexa -1)+1):(6iHexa), 2]
        # 六面体包含的基函数 id
        hexa.inBfsID    .=   @view rbfsID[:,iHexa]
    end
    return rbfsInfo

end

"""
六面体的六个面的端点
"""
const facesVertIDs  =   @SMatrix [  2  1  4  1  5  1;
                                    3  4  3  2  6  2;
                                    7  8  7  6  7  3;
                                    6  5  8  5  8  4]
"""
六面体的六个面的对面
"""
const facen2OppositeID  =   @SVector [2, 1, 4, 3, 6, 5]

"""
六面体的六个面的 对面的 端点
"""                                    
const oppFacesVertIDs   =   facesVertIDs[:, facen2OppositeID]


"""
    getFreeVns(hexa::HexahedraInfo, i::Integer)

得到六面体 `hexa` 第 `i` 个所在的基函数的 所有自由端( ``r₀`` )。
该面定义为 ``uvw`` 坐标中某值为 ``1(0)`` 的面，
自由端应定义为在该面的 ``uvw`` 坐标中，将该值赋值为 ``0(1)`` 时计算得到的点
如 ``u = 1`` 的面为六面体的第 ``(2,6,7,3)`` 个点构成的面，``r₀`` 为 ``uvw`` 
坐标为 ``(0, v, w)`` 的点，此点等同于 ``u = 0`` 的面的四边形的参数坐标为 (v, w) 的点 
在构造 RBF 时面按照 按 ``u=1, u=0, v = 1, v = 0, w = 1, w = 0`` 的顺序排列, 
因此函数根据输入的面的序号进行计算。
"""
function getFreeVns(hexa::HexahedraInfo, i::Integer)
    # 对面的四个点
    @views @inbounds vertices = hexa.vertices[:, oppFacesVertIDs[:, i]]
    vertices * QuadGQInfo.coordinate
end

"""
    getFreeVns(hexa::HexahedraInfo, i::Integer)

得到六面体 `hexa` 处理奇异性时第 `i` 个所在的基函数的 所有自由端( ``r₀`` )。
定义详见[`getFreeVns`](@ref)
"""
function getFreeVnsSglr(hexa::HexahedraInfo, i::Integer)
    # 对面的四个点
    @views @inbounds vertices = hexa.vertices[:, oppFacesVertIDs[:, i]]
    vertices * QuadGQInfoSglr.coordinate
end

"""
    getFreeVnsSSglr(hexa::HexahedraInfo, i::Integer)

得到六面体 `hexa` 处理超奇异性时第 `i` 个所在的基函数的 所有自由端( ``r₀`` )。
定义详见[`getFreeVns`](@ref)
"""
function getFreeVnsSSglr(hexa::HexahedraInfo, i::Integer)
    # 对面的四个点
    @views @inbounds vertices = hexa.vertices[:, oppFacesVertIDs[:, i]]
    vertices * QuadGQInfoSSglr.coordinate
end


"""
    constructFloat2IndexDict(floats::AbstractVector{FT}) where {FT<:AbstractFloat}

构建有序（从小到大）浮点数为键，值为该组数构成的的二重字典，字典值为键的二维数组的线性坐标
"""
function constructFloat2IndexDict(floats::AbstractVector{FT}) where {FT<:AbstractFloat}
    # 预分配结果字典
    re = Dict{FT, Dict{FT, Int64}}()
    # 数组长度
    nF = length(floats)
    # 外层循环
    for ii in 1:nF
        # 内层循环
        relocal = Dict{FT, Int64}()
        for jj in 1:nF
            relocal[floats[jj]] = ii + nF*(jj - 1)
        end # for jj
        re[floats[ii]] = relocal
    end # for ii

    re
end

@doc """
从四边形的 ``uv`` 坐标映射到 自由端的 id
"""
const uv2FreeVnsDict      =   constructFloat2IndexDict(unique(HexaGQInfo.uvw))
const uv2FreeVnsSglrDict  =   constructFloat2IndexDict(unique(HexaGQInfoSglr.uvw))
const uv2FreeVnsSSglrDict =   constructFloat2IndexDict(unique(HexaGQInfoSSglr.uvw))


"""
    selectFreeVnID(uvw::AbstractVector{FT}, i::Integer) where {FT}

根据 ``uvw`` 计算得到第 i 个面所在的基函数的 正常高斯求积 时 “自由端( ``r₀`` )” 的序号。
"""
function selectFreeVnID(uvw::AbstractVector{FT}, i::Integer) where {FT}
    @views @inbounds re = begin
        if (i == 1) || (i == 2)
            uv2FreeVnsDict[uvw[2]][uvw[3]]
        elseif (i == 3) || (i == 4)
            uv2FreeVnsDict[uvw[3]][uvw[1]]
        elseif (i == 5) || (i == 6)
            uv2FreeVnsDict[uvw[1]][uvw[2]]
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    return re
end

"""
    selectFreeVnSglrID(uvw::AbstractVector{FT}, i::Integer) where {FT}

根据 ``uvw`` 得到第 i 个面所在的基函数的 处理奇异性 求积时 “自由端( ``r₀`` )” 的序号。
"""
function selectFreeVnSglrID(uvw::AbstractVector{FT}, i::Integer) where {FT}
    @inbounds re = begin
        if (i == 1) || (i == 2)
            uv2FreeVnsSglrDict[uvw[2]][uvw[3]]
        elseif (i == 3) || (i == 4)
            uv2FreeVnsSglrDict[uvw[3]][uvw[1]]
        elseif (i == 5) || (i == 6)
            uv2FreeVnsSglrDict[uvw[1]][uvw[2]]
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    re
end

"""
    selectFreeVnSSglrID(uvw::AbstractVector{FT}, i::Integer) where {FT}

根据 ``uvw`` 得到第 i 个面所在的基函数的 处理超奇异性 求积时 “自由端( ``r₀`` )” 的序号。
"""
function selectFreeVnSSglrID(uvw::AbstractVector{FT}, i::Integer) where {FT}
    @inbounds re = begin
        if (i == 1) || (i == 2)
            uv2FreeVnsSSglrDict[uvw[2]][uvw[3]]
        elseif (i == 3) || (i == 4)
            uv2FreeVnsSSglrDict[uvw[3]][uvw[1]]
        elseif (i == 5) || (i == 6)
            uv2FreeVnsSSglrDict[uvw[1]][uvw[2]]
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    re
end


"""
    constructGQ1DID2GQ3DIDVector(gqInfo)

构建从六面体体高斯求积点线性索引到三维索引的数组。
"""
function constructGQ1DID2GQ3DIDVector(gqInfo)
    # 求积点坐标
    uvwunique = unique(gqInfo.uvw)
    # 求积点 单个维度上的数量
    gqN1D   =   length(uvwunique)
    # 生成
    return reshape([(i, j, k) for i in 1:gqN1D, j in 1:gqN1D, k in 1:gqN1D], :)
end


"""
构建从六面体的 高斯求积线性索引 映射到 三维索引的数组
"""
const GQ1DID2GQ3DIDVector      =   constructGQ1DID2GQ3DIDVector(HexaGQInfo)
const GQ1DID2GQ3DIDVectorSglr  =   constructGQ1DID2GQ3DIDVector(HexaGQInfoSglr)
const GQ1DID2GQ3DIDVectorSSglr =   constructGQ1DID2GQ3DIDVector(HexaGQInfoSSglr)

"""
    getFreeVIDFromGQ3DID(GQ3DID::NTuple{3, Int}, i::Integer)

得到第 i 个面所在的基函数的 正常高斯求积 下，三维坐标为 `GQ3DID` 的高斯求积点的 “自由端``r₀``” 的序号。
"""
function getFreeVIDFromGQ3DID(GQ3DID::NTuple{3, Int}, i::Integer)
    @inbounds re = begin
        if (i == 1) || (i == 2)
            GQ3DID[2] + (GQ3DID[3] - 1) * GQPNQuad1D
        elseif (i == 3) || (i == 4)
            GQ3DID[1] + (GQ3DID[3] - 1) * GQPNQuad1D
        elseif (i == 5) || (i == 6)
            GQ3DID[1] + (GQ3DID[2] - 1) * GQPNQuad1D
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    re
end

"""
    getFreeVIDFromGQ3DIDSglr(GQ3DID::NTuple{3, Int}, i::Integer)

得到第 i 个面所在的基函数的 处理奇异性时，三维坐标为 `GQ3DID` 的高斯求积点的 “自由端``r₀``” 的序号。
"""
function getFreeVIDFromGQ3DIDSglr(GQ3DID::NTuple{3, Int}, i::Integer)
    @inbounds re = begin
        if (i == 1) || (i == 2)
            GQ3DID[2] + (GQ3DID[3] - 1) * GQPNQuad1DSglr
        elseif (i == 3) || (i == 4)
            GQ3DID[1] + (GQ3DID[3] - 1) * GQPNQuad1DSglr
        elseif (i == 5) || (i == 6)
            GQ3DID[1] + (GQ3DID[2] - 1) * GQPNQuad1DSglr
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    re
end

"""
    getFreeVIDFromGQ3DIDSSglr(GQ3DID::NTuple{3, Int}, i::Integer)

得到第 i 个面所在的基函数的 处理超奇异性时，三维坐标为 `GQ3DID` 的高斯求积点的 “自由端``r₀``” 的序号。
"""
function getFreeVIDFromGQ3DIDSSglr(GQ3DID::NTuple{3, Int}, i::Integer)
    @inbounds re = begin
        if (i == 1) || (i == 2)
            GQ3DID[2] + (GQ3DID[3] - 1) * GQPNQuad1DSSglr
        elseif (i == 3) || (i == 4)
            GQ3DID[1] + (GQ3DID[3] - 1) * GQPNQuad1DSSglr
        elseif (i == 5) || (i == 6)
            GQ3DID[1] + (GQ3DID[2] - 1) * GQPNQuad1DSSglr
        else
            throw("Only 6 faces, $i is a inproper index.")
        end
    end
    re
end