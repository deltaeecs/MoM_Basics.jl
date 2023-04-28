# 本文件构建三角形类和一些积分函数、奇异值处理等
"""TriangleMesh用于保存三角形网格信息：
    trinum: 三角形数量
    node  : 三角形所在节点位置坐标"""
struct TriangleMesh{IT, FT} <: MeshDataType
    trinum      ::Int
    node        ::Array{FT, 2}
    triangles   ::Array{IT, 2}
end

"""
单个三角形信息：\n
triID       ::IT                    编号\n
area        ::FT                    面积\n
verticesID  ::MVector{3, IT}        所在节点id\n
vertices    ::MMatrix{3, 3, FT, 9} 三角形3个角点坐标，每列为一个点\n
center      ::MVec3D{FT}             中心坐标\n
facen̂       ::MVec3D{FT}             面的外法向量\n
edgel       ::MVec3D        三边长\n
edgev̂       ::MMatrix{3, 3, FT, 9} 三个边的指向向量\n
edgen̂       ::MMatrix{3, 3, FT, 9} 三个边的外法向量\n
inBfsID     ::MVector{3, IT}        三角形所在的三个基函数的ID\n
# RWGfreevid  ::MVector{3, IT}     (合理安排位置后，三个基函数的自由端即为三角形三个点的顺序
                                    因此再不需要单独保存)对应三个RWG基函数的自由点在该三角形中的编号\n
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

"""TriangleInfo的默认构造函数，所有元素置零"""
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

# 高斯求积点数
const GQPNTri        =   4
# 处理奇异性时的高斯求积点
const GQPNTriSglr    =   7

~(@isdefined TriGQInfo)         && const TriGQInfo          =   GaussQuadrature4Geo.GaussQuadratureInfo(:Triangle, GQPNTri, Precision.FT)
~(@isdefined TriGQInfoSglr)     && const TriGQInfoSglr      =   GaussQuadrature4Geo.GaussQuadratureInfo(:Triangle, GQPNTriSglr, Precision.FT)


@inline function getGQPTri(tri::TriangleInfo, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tri.vertices * TriGQInfo.coordinate[:, ii]
end

"""
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTri(tri::TriangleInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfo.coordinate
end

@inline function getGQPTriSglr(tri::TriangleInfo, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tri.vertices * TriGQInfoSglr.coordinate[:, ii]
end

""" 
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTriSglr(tri::TriangleInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfoSglr.coordinate
end

## 本想用于创建 zeros函数创建 TriangleInfo 数组， 但结果是对 zeros(TriangleInfo, N) 的结果向量只创建了一个实例，其余全是视图，因此不可用。
# import Base: zero
# zero(::Type{TriangleInfo{IT, FT}}) where{IT<:Integer, FT<:AbstractFloat} = TriangleInfo{IT, FT}()

# TriangleInfo() = TriangleInfo{IntDtype, Precision.FT}()

"""用于在预分配好的数组里写入对应的三角形id和点坐标数据"""
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

"""计算三角形边长、边外法向量、面法向量、面积，直接写入triangleInfo（单个三角形类型实例）"""
function setTriParam!(triangleInfo::TriangleInfo)
    @inbounds @views begin
    # vertices    指向三角形的角点坐标
    vertices            =   triangleInfo.vertices    
    # edgev̂指向三角形的的指向向量(未单位化)
    edgev̂          =   triangleInfo.edgev̂
    # edgev̂[:, 1]    =   vertices[:, 3]  - vertices[:, 2]
    # edgev̂[:, 2]    =   vertices[:, 1]  - vertices[:, 3]
    # edgev̂[:, 3]    =   vertices[:, 2]  - vertices[:, 1]
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

    return

end # function

"""
单个构成四面体的三角形信息：\n
vertices    ::MMatrix{3, 3, FT, 9} 三角形3个角点坐标，每列为一个点\n
edgel       ::MVec3D        三边长\n
edgev̂       ::MMatrix{3, 3, FT, 9} 三个边的指向向量\n
edgen̂       ::MMatrix{3, 3, FT, 9} 三个边的外法向量\n
# RWGfreevid  ::MVector{3, IT}     (合理安排位置后，三个基函数的自由端即为三角形三个点的顺序
                                    因此再不需要单独保存)对应三个RWG基函数的自由点在该三角形中的编号\n
"""
mutable struct Tris4Tetra{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}
    isbd        ::Bool
    δκ          ::Complex{FT}
    vertices    ::MMatrix{3, 3, FT, 9}
    edgel       ::MVec3D{FT}
    edgev̂       ::MMatrix{3, 3, FT, 9}
    edgen̂       ::MMatrix{3, 3, FT, 9}
end

"""TriangleInfo的默认构造函数，所有元素置零"""
function Tris4Tetra{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    isbd        =    true
    δκ          =    zero(Complex{FT})
    vertices    =    zero(MMatrix{3, 3, FT, 9})
    edgel       =    zero(MVec3D{FT})
    edgev̂       =    zero(MMatrix{3, 3, FT, 9})
    edgen̂       =    zero(MMatrix{3, 3, FT, 9})
    return Tris4Tetra{IT, FT}(isbd, δκ, vertices, edgel, edgev̂, edgen̂)
end

"""
计算得到第ii个几何体的高斯求积坐标
    该函数针对三角形， 输入值：
vertices：  三角形角点坐标， MArray{Tuple{3, 3}, FT}
ii      ：  索取第几个坐标
GQMode  :   用于奇异积分时
返回值  :  第ii个高斯求积点坐标, 类型为:SVector{3, FT}
"""
@inline function getGQPTri(tri::Tris4Tetra, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds tri.vertices * TriGQInfo.coordinate[:, ii]
end

"""
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTri(tri::Tris4Tetra)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfo.coordinate
end

@inline function getGQPTriSglr(tri::Tris4Tetra, ii::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds  tri.vertices * TriGQInfoSglr.coordinate[:, ii]
end

""" 
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPTriSglr(tri::Tris4Tetra)
    # 直接向量相乘，此处采用的是矩阵乘法
    @inbounds tri.vertices * TriGQInfoSglr.coordinate
end
