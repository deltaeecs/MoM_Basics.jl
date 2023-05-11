
"""
    Quads4Hexa{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}

单个构成六面体的四边形信息：
```
isbd        ::Bool                  是否在体区域边界，
δκ          ::Complex{FT}           介质对比度变化量，
vertices    ::MMatrix{3, 4, FT, 12} 四边形 4 个角点坐标，每列为一个点，
edgel       ::MVector{4, FT}        四个边长，
edgev̂       ::MMatrix{3, 4, FT, 12} 四个边的单位指向向量，
edgen̂       ::MMatrix{3, 4, FT, 12} 四个边的单位外法向量。
```
合理安排位置后，四个基函数的自由端即为四边形四个点的顺序。
"""
mutable struct Quads4Hexa{IT<: Integer, FT<:AbstractFloat} <: SurfaceCellType{IT, FT}
    isbd        ::Bool
    δκ          ::Complex{FT}
    vertices    ::MMatrix{3, 4, FT, 12}
    edgel       ::MVector{4, FT}
    edgev̂       ::MMatrix{3, 4, FT, 12}
    edgen̂       ::MMatrix{3, 4, FT, 12}
end

"""
    Quads4Hexa{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}

Quads4Hexa的默认构造函数，默认为边界，其它所有元素置零。
"""
function Quads4Hexa{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    isbd        =    true
    δκ          =    zero(Complex{FT})
    vertices    =    zero(MMatrix{3, 4, FT, 12})
    edgel       =    zero(MVector{4, FT})
    edgev̂       =    zero(MMatrix{3, 4, FT, 12})
    edgen̂       =    zero(MMatrix{3, 4, FT, 12})
    return Quads4Hexa{IT, FT}(isbd, δκ, vertices, edgel, edgev̂, edgen̂)
end

"""
    (q::Quads4Hexa)(u::FT, v::FT) where {FT}

计算六面体边界四边形在局部坐标 `(u, v)` 下的点。
"""
function (q::Quads4Hexa)(u::FT, v::FT) where {FT}
    w = SVector{4, FT}((1-u)*(1-v), u*(1-v), u*v, (1-u)*v)
    return q.vertices * w
end

# 改的话和六面体一起改！
"""
六面体边界四边形 1 维度高斯求积点数。
"""
const GQPNQuad1D      =   2
"""
处理奇异性时六面体边界四边形 1 维度高斯求积点数。
"""
const GQPNQuad1DSglr  =   3
"""
处理超奇异性时六面体边界四边形 1 维度高斯求积点数。
"""
const GQPNQuad1DSSglr =   4
"""
六面体边界四边形正常处理高斯求积时高斯求积点数。
"""
const GQPNQuad      =   GQPNQuad1D^2
"""
六面体边界四边形处理奇异性时高斯求积点数。
"""
const GQPNQuadSglr  =   GQPNQuad1DSglr^2
"""
六面体边界四边形处理超奇异性时高斯求积点数。
"""
const GQPNQuadSSglr =   GQPNQuad1DSSglr^2

"""
六面体边界四边形正常处理高斯求积信息。
"""
const QuadGQInfo        =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuad, Precision.FT)
"""
六面体边界四边形处理奇异性高斯求积信息。
"""
const QuadGQInfoSglr    =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuadSglr, Precision.FT)
"""
六面体边界四边形处理超奇异性高斯求积信息。
"""
const QuadGQInfoSSglr   =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuadSSglr, Precision.FT)

"""
    getGQPQuad(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    getGQPQuad(quad::Quads4Hexa)

计算 `quad` 正常求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPQuad(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfo.coordinate[:, ii]
end
function getGQPQuad(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfo.coordinate
end

"""
    getGQPQuadSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    getGQPQuadSglr(quad::Quads4Hexa)

计算 `quad` 处理奇异性的第 `i` 个或所有高斯求积坐标。
"""
function getGQPQuadSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfoSglr.coordinate[:, ii]
end
function getGQPQuadSglr(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfoSglr.coordinate
end

"""
    getGQPQuadSSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    getGQPQuadSSglr(quad::Quads4Hexa)

计算 `quad` 处理超奇异性的第 `i` 个或所有高斯求积坐标。
"""
function getGQPQuadSSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfoSSglr.coordinate[:, ii]
end
function getGQPQuadSSglr(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfoSSglr.coordinate
end