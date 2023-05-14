# 本文件构建四面体形类和一些积分函数、奇异值处理等

"""
    HexahedraMesh{IT, FT} <: MeshDataType

四面体网格文件，包括四面体数 `tetranum`、节点 `node`、构成四面体的节点 id 数组 `tetrahedras` 等信息。
"""
struct TetrahedraMesh{IT, FT} <: MeshDataType
    tetranum      ::Int
    node          ::Array{FT, 2}
    tetrahedras   ::Array{IT, 2}
end

"""
    TetrahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex} <: VolumeCellType{IT, FT, CT}

四面体信息结构体：
```
tetraID     ::IT                    编号
volume      ::FT                    体积
ε           ::CT                    相对介电常数
κ           ::CT                    介质对比度
center      ::MVec3D{FT}            中心坐标
verticesID  ::MVector{4, IT}        所在节点id
vertices    ::MMatrix{3, 4, FT, 12} 四面体4个角点坐标，每列为一个点
facesn̂      ::MMatrix{3, 4, FT, 12} 四个面的外法向量
facesvid    ::MMatrix{3, 4, IT, 12} 四个面包含的三个id
facesArea   ::MVector{4, FT}        四个面的面积（根据为unitri正负部分赋予正负号）
faces       ::Vector{Tris4Tetra{IT, FT}}    四个面的具体信息
inBfsID     ::Vector{IT}            四面体所在的基函数的ID
```
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
    TetrahedraInfo{IT, FT, CT}(hexaID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat, CT<:Complex}

`TetrahedraInfo` 的默认构造函数，除了输入的编号 `tetraID` 外所有元素置零。
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


```
正常求积四面体高斯求积点数。
```
const GQPNTetra      =   5
```
处理奇异性时的四面体高斯求积点数。
```
const GQPNTetraSglr  =   11

```
正常求积四面体高斯求积信息。
```
const TetraGQInfo        =   GaussQuadrature4Geo.GaussQuadratureInfo(:Tetrahedron, GQPNTetra, Float32)
```
处理奇异性时四面体高斯求积信息。
```
const TetraGQInfoSglr    =   GaussQuadrature4Geo.GaussQuadratureInfo(:Tetrahedron, GQPNTetraSglr, Float32)

"""
    setHexaCoor!( tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}, tetraMeshData::TetrahedraMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

在预分配好的四面体数组 `tetrasInfo` 里写入 `tetraMeshData` 中对应的四面体编号、点坐标、中心位置数据。
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
    setTetraParam!(tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

计算四面体体积、面外法向量、面积，并写入 `tetrasInfo` 。
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
    getTetrasInfo(tetrameshData::TetrahedraMesh{IT, FT}, VolumeBFType::Symbol) where{IT, FT}


根据四面体网格信息 `tetrameshData` 和体基函数类型 `VolumeBFType` 生成网格信息向量 `tetrasInfo` 和基函数信息向量 `bfsInfo` 。
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
    setδκ!(geosInfo::AbstractVector{VT}) where {VT<:VolumeCellType}

设置体网格信息 `geosInfo` 中每个面上的介质对比度差值。
"""
function setδκ!(geosInfo::AbstractVector{VT}) where {VT<:VolumeCellType}
    # 循环置零
    @threads for it in eachindex(geosInfo)
        # 体网格
        geo   =   geosInfo[it]
        # 对面循环
        for face in geo.faces
            # 置零
            face.δκ = 0
        end
    end
    # 循环设置δκ
    @threads for it in eachindex(geosInfo)
        # 体网格
        geo   =   geosInfo[it]
        # 四面体的介质对比度
        κ   =    geo.κ
        # 对面循环
        for iface in eachindex(geo.faces)
            # 面
            face =  geo.faces[iface]
            # 根据在正负面决定加上或减去 κ
            temp =  (geo.facesArea[iface] > 0 ? κ : -κ)
            # 写入值
            face.δκ += temp
        end
    end
end


@doc """
    getGQPTetra(tetra::TetrahedraInfo, i::IT) where {IT <: Integer}
    getGQPTetra(tetra::TetrahedraInfo)

计算 `tetra` 正常求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTetra(tetra::TetrahedraInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tetra.vertices * TetraGQInfo.coordinate[:, i]
end
function getGQPTetra(tetra::TetrahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tetra.vertices * TetraGQInfo.coordinate
end

@doc """
    getGQPTetraSglr(tetra::TetrahedraInfo, i::IT) where {IT <: Integer}
    getGQPTetraSglr(tetra::TetrahedraInfo)

计算 `tetra` 处理奇异性求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTetraSglr(tetra::TetrahedraInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tetra.vertices * TetraGQInfoSglr.coordinate[:, i]
end
function getGQPTetraSglr(tetra::TetrahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tetra.vertices * TetraGQInfoSglr.coordinate
end