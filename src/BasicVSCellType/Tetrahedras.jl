# 本文件构建四面体形类和一些积分函数、奇异值处理等

"""
四面体网格文件，包括四面体数、节点、构成四面体的节点 id 数组等信息
"""
struct TetrahedraMesh{IT, FT} <: MeshDataType
    tetranum      ::Int
    node          ::Array{FT, 2}
    tetrahedras   ::Array{IT, 2}
end


"""
单个四面体信息：\n
tetraID     ::IT                    编号\n
volume      ::FT                    体积\n
ε           ::CT                    相对介电常数\n
κ           ::CT                    介质对比度\n
center      ::MVec3D{FT}            中心坐标\n
verticesID  ::MVector{4, IT}        所在节点id\n
vertices    ::MMatrix{3, 4, FT, 12} 四面体4个角点坐标，每列为一个点\n
facesn̂      ::MMatrix{3, 4, FT, 12} 四个面的外法向量\n
facesvid    ::MMatrix{3, 4, IT, 12} 四个面包含的三个id
facesArea   ::MVector{4, FT}        四个面的面积（根据为unitri正负部分赋予正负号）\n
faces       ::Vector{Tris4Tetra{IT, FT}}    四个面的具体信息\n
inBfsID     ::Vector{IT}            四面体所在的基函数的ID\n
"""
mutable struct TetrahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex} <: VolumeCellType{IT, FT, CT}
    tetraID     ::IT
    volume      ::FT
    ε           ::CT
    κ           ::CT
    center      ::MVec3D{FT}
    verticesID  ::MVector{4, IT}
    vertices    ::MMatrix{3, 4, FT, 12}
    facesn̂      ::MMatrix{3, 4, FT, 12}
    facesArea   ::MVector{4, FT}
    faces       ::Vector{Tris4Tetra{IT, FT}}
    inBfsID     ::Vector{IT}
end

"""
TriangleInfo的默认构造函数，所有元素置零\n
"""
function TetrahedraInfo{IT, FT, CT}(tetraID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat, CT<:Complex}
    volume::FT  =    zero(FT)
    ε           =    one(CT)
    κ           =    zero(CT)
    center      =    zero(MVec3D{FT})
    verticesID  =    zero(MVector{4, IT})
    vertices    =    zero(MMatrix{3, 4, FT, 12})
    facesn̂      =    zero(MMatrix{3, 4, FT, 12})
    facesArea   =    zero(MVector{4, FT})
    faces       =    Vector{Tris4Tetra{IT, FT}}(undef, 4)
    inBfsID     =    zeros(IT, 4)
    return TetrahedraInfo{IT, FT, CT}(tetraID,  volume, ε, κ, center, verticesID, vertices,
                                        facesn̂, facesArea, faces, inBfsID)
end


# 高斯求积点数
const GQPNTetra      =   5
# 处理奇异性时的高斯求积点
const GQPNTetraSglr  =   11


~(@isdefined TetraGQInfo)       && const TetraGQInfo        =   GaussQuadrature4Geo.GaussQuadratureInfo(:Tetrahedron, GQPNTetra, Precision.FT)
~(@isdefined TetraGQInfoSglr)   && const TetraGQInfoSglr    =   GaussQuadrature4Geo.GaussQuadratureInfo(:Tetrahedron, GQPNTetraSglr, Precision.FT)

"""
用于在预分配好的数组里写入对应的四面体id和点坐标数据
"""
function setTetraCoor!( tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}, 
                        tetraMeshData::TetrahedraMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}
    @threads for i in 1:length(tetrasInfo)
        tetrasInfo[i].tetraID      =  i
        tetrasInfo[i].verticesID  .=  view(tetraMeshData.tetrahedras, :, i)
        for ii in 1:4
            tetrasInfo[i].vertices[:,ii]   .=  @view tetraMeshData.node[:,tetrasInfo[i].verticesID[ii]]
        end
        tetrasInfo[i].center   =   mean(tetrasInfo[i].vertices, dims = 2)
    end

    return nothing
end

"""
计算四面体边长、面外法向量、面积，直接写入tetraInfo（单个四面体类型实例）
"""
function setTetraParam!(tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

    vertloop    =   Int[1, 2, 3, 4, 1, 2, 3]
    @threads for i in 1:length(tetrasInfo)
        @inbounds begin
            # 第 i 四面体
            tetraInfo    =   tetrasInfo[i]
            # 四面体四个点的坐标
            vertices    =   tetraInfo.vertices

            # 计算四面体的某个点（第一个）指向第其他个点的三个向量
            @views vert12  =    vertices[:, 2] .- vertices[:, 1]
            @views vert13  =    vertices[:, 3] .- vertices[:, 1]
            @views vert14  =    vertices[:, 4] .- vertices[:, 1]
            # 四面体体积
            tetraInfo.volume    =   cross(vert12, vert13) ⋅ vert14 / 6
            
            # 对四个面循环计算面的外法向量，注意第 facei 个面为 第 facei 个 vertice 对面的三个点构成的四面体
            for facei in 1:4
                #  构成面的三个点的id verts ids
                faceiIDs    =    view(vertloop, (facei + 1):(facei + 3))
                # 三个点的视图
                faceVertViews   =   view(vertices, :, faceiIDs)
                # 计算面法向量（未定向\未归一化）
                @views facen̂    =   tetraInfo.facesn̂[:, facei]
                @views facen̂   .=   cross(faceVertViews[:, 2] .- faceVertViews[:, 1],  faceVertViews[:, 3] .- faceVertViews[:, 1])
                # 未归一化面法向量的模即为四面体面积两倍
                faceArea    =   norm(facen̂)/2
                # 归一化面法向量
                facen̂     ./=   2faceArea

                # 定向，第 facei 个 vertice 点到其余三个点中的一个的向量与外法向量的点乘结果应该大于0
                @views ((faceVertViews[:, 1] .- vertices[:, facei]) ⋅ facen̂ < 0) && begin facen̂ .*= -1 end

                # 写入结果
                tetraInfo.facesArea[facei]  =   faceArea
                # 面法向量结果采用的是内存视图，不必再重新写入
            end

        end
    end
    return

end # function


"""
将四面体、基函数相关的构造函数封装在此函数里
"""
function getTetrasInfo(tetrameshData::TetrahedraMesh{IT, FT}, VolumeBFType::Symbol) where{IT, FT}
    # 分配四面体信息类型数组
    tetrasInfo       =   [TetrahedraInfo{IT, FT, Complex{FT}}()  for _ in 1:tetrameshData.tetranum ]
    # 写入对应的四面体id和点坐标数据
    setTetraCoor!(tetrasInfo, tetrameshData)
    # 写入四面体边长、边外法向量、面法向量、面积
    setTetraParam!(tetrasInfo)
    # 计算基函数信息并修改四面体信息中基函数相关项
    bfsInfo  =   setTriangles4Tetras!(tetrameshData, tetrasInfo, Val(VolumeBFType))
    return tetrasInfo, bfsInfo
end

"""
此函数用于设置 SWG 基函数的两个四面体的介质对比度差值
"""
function setδκ!(tetrasInfo::AbstractVector{TetrahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 线程锁，防止对同一数据操作引起的冲突
    κLock = SpinLock()
    # 循环置零
    @threads for it in eachindex(tetrasInfo)
        # 四面体
        tetra   =   tetrasInfo[it]
        # 对面循环
        for iface in 1:4
            # 面
            face =  tetra.faces[iface]
            # 置零
            face.δκ = 0
        end
    end
    # 循环设置δκ
    @threads for it in eachindex(tetrasInfo)
        # 四面体
        tetra   =   tetrasInfo[it]
        # 四面体的介质对比度
        κ   =    tetra.κ
        # 对面循环
        for iface in 1:4
            # 面
            face =  tetra.faces[iface]
            # 根据在正负面决定加上或减去 κ
            temp =  (tetra.facesArea[iface] > 0 ? κ : -κ)
            # 写入值
            lock(κLock)
            face.δκ += temp
            unlock(κLock)
        end
    end
end


"""
计算得到第ii个几何体的高斯求积坐标
    该函数针对四面体， 输入值：
vertices：  四面体角点坐标， MArray{Tuple{3, 3}, FT}
ii      ：  索取第几个坐标
GQMode  :   用于奇异积分时
返回值  :  第ii个高斯求积点坐标, 类型为:SVector{3, FT}
"""
@inline function getGQPTetra(tetra::TetrahedraInfo, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tetra.vertices * TetraGQInfo.coordinate[:, ii]
end

"""
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTetra(tetra::TetrahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tetra.vertices * TetraGQInfo.coordinate
end

@inline function getGQPTetraSglr(tetra::TetrahedraInfo, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tetra.vertices * TetraGQInfoSglr.coordinate[:, ii]
end

"""
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTetraSglr(tetra::TetrahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tetra.vertices * TetraGQInfoSglr.coordinate
end