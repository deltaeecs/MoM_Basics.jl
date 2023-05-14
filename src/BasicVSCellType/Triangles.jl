# 本文件构建三角形类和一些积分函数、奇异值处理等
"""
    TriangleMesh{IT, FT} <: MeshDataType

三角形网格文件，包括三角形数 `trinum`、节点 `node`、构成三角形的节点 id 数组 `triangles` 等信息。
"""
struct TriangleMesh{IT, FT} <: MeshDataType
    trinum      ::Int
    node        ::Array{FT, 2}
    triangles   ::Array{IT, 2}
end

"""
    TriangleInfo{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}

三角形信息结构体：
```
triID       ::IT                    编号
area        ::FT                    面积
verticesID  ::MVector{3, IT}        所在节点id
vertices    ::MMatrix{3, 3, FT, 9}  三角形3个角点坐标，每列为一个点
center      ::MVec3D{FT}            中心坐标
facen̂       ::MVec3D{FT}            面的外法向量
edgel       ::MVec3D                三边长
edgev̂       ::MMatrix{3, 3, FT, 9}  三个边的指向向量
edgen̂       ::MMatrix{3, 3, FT, 9}  三个边的外法向量
inBfsID     ::MVector{3, IT}        三角形所在的三个基函数的ID
```
合理安排位置后，三个基函数的自由端即为三角形三个点的顺序。
"""
mutable struct TriangleInfo{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}
    triID       ::IT
    ε           ::Complex{FT}
    area        ::FT
    verticesID  ::MVec3D{IT}
    vertices    ::MMatrix{3, 3, FT, 9}
    center      ::MVec3D{FT}
    facen̂       ::MVec3D{FT}
    edgel       ::MVec3D{FT}
    edgev̂       ::MMatrix{3, 3, FT, 9}
    edgen̂       ::MMatrix{3, 3, FT, 9}
    inBfsID     ::MVec3D{IT}
end

"""
    TriangleInfo{IT, FT}(triID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat}

`TriangleInfo` 的默认构造函数，除了输入的编号 `triID` 外所有元素置零。
"""
function TriangleInfo{IT, FT}(triID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat}
    ε           =    one(Complex{FT})*ε_0
    area::FT    =    zero(FT)
    verticesID  =    zero(MVec3D{IT})
    vertices    =    zero(MMatrix{3, 3, FT, 9})
    center      =    zero(MVec3D{FT})
    facen̂       =    zero(MVec3D{FT})
    edgel       =    zero(MVec3D{FT})
    edgev̂       =    zero(MMatrix{3, 3, FT, 9})
    edgen̂       =    zero(MMatrix{3, 3, FT, 9})
    inBfsID     =    zero(MVec3D{IT})
    return TriangleInfo{IT, FT}(triID,  ε,  area,  verticesID, vertices    ,
                                center, facen̂, edgel, edgev̂, 
                                edgen̂, inBfsID)
end

```
正常求积三角形高斯求积点数。
```
const GQPNTri        =   4
```
处理奇异性时的三角形高斯求积点数。
```
const GQPNTriSglr    =   7

```
正常求积三角形高斯求积信息。
```
const TriGQInfo          =   GaussQuadrature4Geo.GaussQuadratureInfo(:Triangle, GQPNTri, Float32)
```
处理奇异性时三角形高斯求积信息。
```
const TriGQInfoSglr      =   GaussQuadrature4Geo.GaussQuadratureInfo(:Triangle, GQPNTriSglr, Float32)

@doc """
    getGQPTri(tri::TriangleInfo, i::IT) where {IT <: Integer}
    getGQPTri(tri::TriangleInfo)

计算 `tri` 正常求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTri(tri::TriangleInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tri.vertices * TriGQInfo.coordinate[:, i]
end
function getGQPTri(tri::TriangleInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfo.coordinate
end
@doc """
    getGQPTriSglr(tri::TriangleInfo, i::IT) where {IT <: Integer}
    getGQPTriSglr(tri::TriangleInfo)

计算 `tri` 处理奇异性求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTriSglr(tri::TriangleInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tri.vertices * TriGQInfoSglr.coordinate[:, i]
end
@inline function getGQPTriSglr(tri::TriangleInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfoSglr.coordinate
end

"""
    setTricoor!( trianglesInfo::Vector{TriangleInfo{IT, FT}}, TriangleMeshData::TriangleMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat}

在预分配好的三角形数组 `trianglesInfo` 里写入 `TriangleMeshData` 中对应的三角形编号、点坐标、中心位置数据。
"""
function setTricoor!(   trianglesInfo::Vector{TriangleInfo{IT, FT}}, 
                        TriangleMeshData::TriangleMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat}
    @views @inbounds for i in 1:length(trianglesInfo)
        trianglesInfo[i].triID           =  i
        trianglesInfo[i].verticesID     .=  TriangleMeshData.triangles[:, i]
        for ii in 1:3
            trianglesInfo[i].vertices[:,ii]   =  TriangleMeshData.node[:,trianglesInfo[i].verticesID[ii]]
        end
    end

    return nothing
end

"""
    setTriParam!(triangleInfo::TriangleInfo)

计算三角形边长、边外法向量、面法向量、面积，直接写入 `triangleInfo` 。
"""
function setTriParam!(triangleInfo::TriangleInfo)
    @inbounds @views begin
    # vertices    指向三角形的角点坐标
    vertices            =   triangleInfo.vertices    
    # edgev̂指向三角形的的指向向量(未单位化)
    edgev̂          =   triangleInfo.edgev̂
    @views for edgei in 1:3
        edgev̂[:, edgei]    =   vertices[:, EDGEVpINTriVsID[edgei]]  - vertices[:, EDGEVmINTriVsID[edgei]]
    end

    # 边长
    for edgei in 1:3
        @views triangleInfo.edgel[edgei]   =   norm(edgev̂[:, edgei])
    end

    # 面法向量（未单位化）
    @views triangleInfo.facen̂[:]   =  -cross(edgev̂[:, 1], edgev̂[:, 3])

    # 面积
    triangleInfo.area       =   norm(triangleInfo.facen̂)/2

    # 中心
    triangleInfo.center     =   mean(vertices, dims=2)

    # 面法向量单位化
    triangleInfo.facen̂    /=   2*triangleInfo.area

    
    # 计算边的指向向量、外法向量
    edgen̂                  =   triangleInfo.edgen̂
    for edgei in 1:3
        # 边的指向向量单位化
        @views edgev̂[:, edgei]   /=   triangleInfo.edgel[edgei]
        # 最后计算边外法向量
        @views edgen̂[:, edgei]    =   cross(edgev̂[:, edgei], triangleInfo.facen̂)
    end

    end # begin

    return nothing

end # function

"""
    TriangleInfo{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}

构成四面体的三角形信息结构体：
```
isbd        ::Bool                  是否在边界上
δκ          ::Complex{FT}           面两侧介质对比度差值
vertices    ::MMatrix{3, 3, FT, 9}  三角形3个角点坐标，每列为一个点
edgel       ::MVec3D{FT}            三边长
edgev̂       ::MMatrix{3, 3, FT, 9}  三个边的指向向量
edgen̂       ::MMatrix{3, 3, FT, 9}  三个边的外法向量
```
合理安排位置后，三个基函数的自由端即为三角形三个点的顺序。
"""
mutable struct Tris4Tetra{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}
    isbd        ::Bool
    δκ          ::Complex{FT}
    vertices    ::MMatrix{3, 3, FT, 9}
    edgel       ::MVec3D{FT}
    edgev̂       ::MMatrix{3, 3, FT, 9}
    edgen̂       ::MMatrix{3, 3, FT, 9}
end

"""
    Tris4Tetra{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}

`Tris4Tetra` 的默认构造函数，默认在边界上，其它所有元素置零。
"""
function Tris4Tetra{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    isbd        =    true
    δκ          =    zero(Complex{FT})
    vertices    =    zero(MMatrix{3, 3, FT, 9})
    edgel       =    zero(MVec3D{FT})
    edgev̂       =    zero(MMatrix{3, 3, FT, 9})
    edgen̂       =    zero(MMatrix{3, 3, FT, 9})
    return Tris4Tetra{IT, FT}(isbd, δκ, vertices, edgel, edgev̂, edgen̂)
end

@doc """
    getGQPTri(tri::Tris4Tetra, i::IT) where {IT <: Integer}
    getGQPTri(tri::Tris4Tetra)

计算 `tri` 正常求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTri(tri::Tris4Tetra, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tri.vertices * TriGQInfo.coordinate[:, i]
end
function getGQPTri(tri::Tris4Tetra)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfo.coordinate
end

@doc """
    getGQPTriSglr(tri::Tris4Tetra, i::IT) where {IT <: Integer}
    getGQPTriSglr(tri::Tris4Tetra)

计算 `tri` 处理奇异性求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPTriSglr(tri::Tris4Tetra, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tri.vertices * TriGQInfoSglr.coordinate[:, i]
end
function getGQPTriSglr(tri::Tris4Tetra)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfoSglr.coordinate
end
