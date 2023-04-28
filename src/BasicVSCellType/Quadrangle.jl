
"""
单个构成六面体的四边形信息：\n
vertices    ::MMatrix{3, 4, FT, 12} 四边形 4 个角点坐标，每列为一个点\n
edgel       ::MVector{4, FT}       四边长\n
edgev̂       ::MMatrix{3, 4, FT, 12} 四个边的指向向量\n
edgen̂       ::MMatrix{3, 4, FT, 12} 四个边的外法向量\n
# RWGfreevid  ::MVector{3, IT}     (合理安排位置后，四个基函数的自由端即为四边形四个点的顺序
                                    因此再不需要单独保存)对应四个RWG基函数的自由点在该四边形中的编号\n
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
Quads4Hexa的默认构造函数，所有元素置零
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
Quads4Hexa的取高斯求积点函数
"""
function (q::Quads4Hexa)(u::FT, v::FT) where {FT}
    w = SVector{4, FT}((1-u)*(1-v), u*(1-v), u*v, (1-u)*v)
    return q.vertices * w
end

# 改的话和六面体一起改！
# 1维度高斯求积点数
const GQPNQuad1D      =   2
# 处理奇异性时的高斯求积点
const GQPNQuad1DSglr  =   3
# 处理超奇异性时的高斯求积点
const GQPNQuad1DSSglr =   4
# 高斯求积点数
const GQPNQuad      =   GQPNQuad1D^2
# 处理奇异性时的高斯求积点
const GQPNQuadSglr  =   GQPNQuad1DSglr^2
# 处理超奇异性时的高斯求积点
const GQPNQuadSSglr =   GQPNQuad1DSSglr^2


~(@isdefined QuadGQInfo)       && const QuadGQInfo        =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuad, Precision.FT)
~(@isdefined QuadGQInfoSglr)   && const QuadGQInfoSglr    =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuadSglr, Precision.FT)
~(@isdefined QuadGQInfoSSglr)  && const QuadGQInfoSSglr   =   GaussQuadrature4Geo.GaussQuadratureInfo(:Quadrangle, GQPNQuadSSglr, Precision.FT)


"""
计算得到第ii个几何体的高斯求积坐标
    该函数针对四边形， 输入值：
vertices：  四边形角点坐标， MArray{Tuple{3, 3}, FT}
ii      ：  索取第几个坐标
GQMode  :   用于奇异积分时
返回值  :  第ii个高斯求积点坐标, 类型为:SVector{3, FT}
"""
@inline function getGQPQuad(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfo.coordinate[:, ii]
end

"""
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPQuad(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfo.coordinate
end

@inline function getGQPQuadSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfoSglr.coordinate[:, ii]
end

""" 
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPQuadSglr(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfoSglr.coordinate
end

@inline function getGQPQuadSSglr(quad::Quads4Hexa, ii::IT) where {IT <: Integer}
    @views @inbounds quad.vertices * QuadGQInfoSSglr.coordinate[:, ii]
end

""" 
同上函数，此时不输入ii, 返回所有求积点，返回值类型为
"""
@inline function getGQPQuadSSglr(quad::Quads4Hexa)
    @views @inbounds quad.vertices * QuadGQInfoSSglr.coordinate
end