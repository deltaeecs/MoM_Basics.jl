# 定义SWG基函数数据结构
"""
SWG 基函数复合类型参数解释：
```
isbd        ::Bool              是否为边界元即半基函数
bfID        ::IT                基函数编号，整形
inGeo       ::MVector{2, IT}    基函数所在两个四面体编号（半基函数为1个，赋值成一样的两个），长度为2的向量数组
inGeoID     ::MVector{2, IT}    基函数在两个四面体中的局部编号（半基函数为1个，赋值成一样的两个），取值1:4，长度为2的向量数组
center      ::MVec3D{FT}        基函数公共面中心，用于八叉树分组
```
"""
mutable struct SWG{IT<:Integer , FT<:AbstractFloat} <: LinearBasisFunction
    isbd        ::Bool
    bfID        ::IT
    inGeo       ::MVector{2, IT}
    inGeoID     ::MVector{2, IT}
    center      ::MVec3D{FT}
end

"""
    SWG{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}

SWG 的默认构造函数，所有元素置零。
"""
function SWG{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    isbd      =    true
    bfID      =    zero(IT)
    inGeo     =    zero(MVector{2, IT})
    inGeoID   =    zero(MVector{2, IT})
    center    =    zero(MVec3D{FT})
    return SWG{IT, FT}(isbd, bfID, inGeo, inGeoID, center)
end
SWG()   =    SWG{IntDtype, Precision.FT}()


"""
    setTriangles4Tetras!(tetrameshData::TetrahedraMesh{IT, FT}, tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}, ::Val{:SWG}) where {IT, FT, CT}

计算构成四面体的所有三角形，并将这些信息写入四面体 `tetrasInfo`，给 SWG 基函数赋值。
"""
function setTriangles4Tetras!(tetrameshData::TetrahedraMesh{IT, FT}, tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}, ::Val{:SWG}) where {IT, FT, CT}
    # 四面体数
    tetranum  =   tetrameshData.tetranum
    # 创建数组，用于后续保存所有面的点和面的对点信息(4)，然后是面所在的四面体id(1)，面所在 unitri(独立三角形) 编号(1)、正负等信息(2)，共 8 列
    face3opp1tri2unitri2   =   zeros(IT, (4tetrameshData.tetranum, 8))
    # 利用内存视图提取其中的所有面的点和面的对点信息信息
    face3opp1           =   reshape(view(face3opp1tri2unitri2, :, 1:4), (:, 16))
    # 利用内存视图提取四面体信息的tri2和基函数信息的unitri2
    tetra2              =   @view face3opp1tri2unitri2[:, 5:6]
    unitri2             =   @view face3opp1tri2unitri2[:, 7:8]

    ## 创建四面体索引( id+localid )并根据排序更新
    # 此维度为四面体id
    tetra2[:, 1]  .=   repeat(1:tetranum,  outer = 4 )
    # 此维度为该面在所在四面体的局部id
    tetra2[:, 2]  .=   repeat(1:4,  inner = tetranum )
    
    # 将每个四面体的四个面的点及其对点存放在 face3opp1 数组中
    #* 注意不可改变下面的顺序因为跟RWG基函数对应的自由点在三角形中顺序有关 *#
    face3opp1[:, 1] .=   face3opp1[:, 8] .=   face3opp1[:, 11] .=   face3opp1[:, 14]  .=  tetrameshData.tetrahedras[2, :]
    face3opp1[:, 2] .=   face3opp1[:, 5] .=   face3opp1[:, 12] .=   face3opp1[:, 15]  .=  tetrameshData.tetrahedras[3, :]
    face3opp1[:, 3] .=   face3opp1[:, 6] .=   face3opp1[:,  9] .=   face3opp1[:, 16]  .=  tetrameshData.tetrahedras[4, :]
    face3opp1[:, 4] .=   face3opp1[:, 7] .=   face3opp1[:, 10] .=   face3opp1[:, 13]  .=  tetrameshData.tetrahedras[1, :]
    # 重塑形状将每个边点组合放入同一维度
    face3opp1        =   reshape(face3opp1, (:, 4))

    ##  对前三列排序使得同一面的表示形式相同
    #   对前三个构成边的点排序
    @views sort!(face3opp1tri2unitri2[:, 1:3], dims = 2)
    # 将每个四面体的四个角标点作为整体排序，次序按v1,v2,v3
    copyto!(face3opp1tri2unitri2, sortslices(face3opp1tri2unitri2, dims = 1, by = x -> (x[1], x[2], x[3]), alg = ThreadsX.MergeSort))

    ## for 循环计算RWG基函数id
    @inbounds let unitriplusid   =   1
        unitri2[1, 1]  =   1; unitri2[1, 2] = 1
        # 初始值即1的位置要预计算一下
        unitriplusid   +=  1

        # 从2开始循环，跟前一边不同为 unitri+，跟前一边不同为 unitri-
        for ii in 2:(tetranum*4)
            if  (face3opp1[ii, 1]  != face3opp1[ii-1, 1]) ||
                (face3opp1[ii, 2]  != face3opp1[ii-1, 2]) ||
                (face3opp1[ii, 3]  != face3opp1[ii-1, 3])
                unitri2[ii, 1]      =   unitriplusid
                unitri2[ii, 2]      =   1
                unitriplusid       +=   1
            elseif  (face3opp1[ii, 1]  ==   face3opp1[ii-1, 1]) &&
                (face3opp1[ii, 2]  ==   face3opp1[ii-1, 2]) && 
                (face3opp1[ii, 3]  ==   face3opp1[ii-1, 3])
                unitri2[ii, 1]      =   unitriplusid  - 1
                unitri2[ii, 2]      =   -1
            end #if
        end #for
    end # let

    # unitri函数数量( 同时也是独立三角形的数量 )
    unitrinum   =   unitri2[end, 1]
    # 单独三角形所在的 face3opp1 的行 id
    unitri2nodeRow    =   zeros(IT, unitrinum)
    iunitri     =   1
    for inode in 1:(tetranum*4)
        (unitri2[inode, 2] == 1) && begin
            unitri2nodeRow[iunitri]     =   inode
            iunitri += 1
        end
    end
    # 预分配所有的构成四面体的三角形
    uniqueTris  =   [Tris4Tetra{IT, FT}()  for _ in 1:unitrinum ]

    # 设置 uniqueTris 的三个点、边长、边法向量等信息
    @threads for iTri in 1:unitrinum
        # 该三角形
        tri =   uniqueTris[iTri]
        vertices    = tri.vertices

        # 该三角形对应的三个角点id
        trivertIDs  =   view(face3opp1, unitri2nodeRow[iTri], :)
        # 设置节点
        for iv in 1:3
            vertices[:, iv]     .=   view(tetrameshData.node, :, trivertIDs[iv])
        end

        # 边指向向量
        edgev̂          =   tri.edgev̂
        @views for edgei in 1:3
            edgev̂[:, edgei]    =   vertices[:, EDGEVpINTriVsID[edgei]]  - vertices[:, EDGEVmINTriVsID[edgei]]
        end

        # 边长
        for edgei in 1:3
            @views tri.edgel[edgei]   =   norm(edgev̂[:, edgei])
        end

        @views facen̂    =  -cross(edgev̂[:, 1], edgev̂[:, 3])
        facen̂    ./=   norm(facen̂)
        # 计算边的指向单位向量、外法单位向量
        edgen̂                  =   tri.edgen̂
        for edgei in 1:3
            # 边的指向向量单位化
            @views edgev̂[:, edgei]  ./=   tri.edgel[edgei]
            # 最后计算边外法向量
            @views edgen̂[:, edgei]    =   cross(edgev̂[:, edgei], facen̂)
        end

    end

    ## 开始构造 SWG 基函数
    # swg函数数量
    swg2        =   unitri2
    swgnum      =   swg2[end, 1]
    # 保存RWG基函数信息的数组
    swgsInfo    =   [SWG{IT, FT}() for _ in 1:swgnum]
    # 写入公共面数是否为边界元信息、基函数相关信息
    @threads for faceii in 1:4tetranum
        # 面id
        faceID      =   swg2[faceii, 1]
        # 正负？
        ispOrm      =   swg2[faceii, 2]
        # 对应基函数
        swgInfo     =   swgsInfo[faceID]
        # 对应面（四边形）
        uniqueTri   =   uniqueTris[faceID]
        # 基函数序号
        swgInfo.bfID    =   faceID
        # 基函数中心
        swgInfo.center[:]   =   mean(uniqueTri.vertices, dims = 2)
        # 其它信息
        (ispOrm == 1) && begin
            #  + 基函数所在六面体、是该六面体第几个基函数
            swgInfo.inGeo[1]        =   tetra2[faceii,1]
            swgInfo.inGeoID[1]      =   tetra2[faceii,2]
        end
        (ispOrm == -1) && begin
            # 有 - 基函数的必为非边界元
            uniqueTri.isbd          =   false
            swgInfo.isbd            =   false
            #  - 基函数所在六面体、是该六面体第几个基函数
            swgInfo.inGeo[2]        =   tetra2[faceii,1]
            swgInfo.inGeoID[2]      =   tetra2[faceii,2]
        end
    end

    ## 设置四面体所包含的基函数信息
    # 按四面体编号排列，效果是相邻的四个三角形同属一个四面体并桉顺序排列
    copyto!(face3opp1tri2unitri2, sortslices(face3opp1tri2unitri2, dims = 1, by = x -> (x[5], x[6]), alg = ThreadsX.MergeSort))
    swgsID  =   reshape(swg2[:, 1], (4,:))
    # 找出每个四面体对应的四个面，并在 “-” 面将三角形面积变负
    @threads for iTera in eachindex(tetrasInfo)

        # 该四面体
        tetra = tetrasInfo[iTera]
        # 四个三角形 id 
        tetra.faces      .=   uniqueTris[unitri2[(4(iTera -1)+1):(4iTera), 1]]
        # 将在负半部分的面积取反
        tetra.facesArea .*=   @view unitri2[(4(iTera -1)+1):(4iTera), 2]
        # 四面体包含的基函数 id
        tetra.inBfsID     =   @view swgsID[:,iTera]

    end
    # tetraDataUnion(tetrasInfo, uniqueTris)
    return swgsInfo

end