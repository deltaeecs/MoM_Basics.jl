"""创建模组用于支持集合体高斯求积运算。"""
module GaussQuadrature4Geo

using FastGaussQuadrature

export GaussQuadratureInfo

using StaticArrays


"""创建类型用于保存高斯求积相对坐标、权重等信息:
    类型参数：
    GeoN: 3代表三角形，4代表四面体
    GQN : 高斯求积点数"""
struct GaussQuadratureInfoStruct{FT<:Real}
    coordinate  ::SMatrix
    weight      ::SVector
    uvw         ::SMatrix
end

"""
GaussQuadratureInfo
构造函数
imput
GeoS :: :Triangle, :Quadrangle, :Tetrahedron, :Hexahedron
"""
function GaussQuadratureInfo(GeoS::Symbol, GQN::IT, FT::DataType = Float64) where {IT<:Integer}
    uvw = SMatrix{0, 0, FT, 0}()
    # 三角形
    if GeoS == :Triangle
        coordinate, weight = gaussQuadratureTri(GQN, FT)
        return GaussQuadratureInfoStruct{FT}(coordinate, weight, uvw)
    # 四面体
    elseif GeoS == :Quadrangle
        coordinate, weight, uvw = gaussQuadratureQuad(GQN, FT)
        return GaussQuadratureInfoStruct{FT}(coordinate, weight, uvw)
    elseif GeoS == :Tetrahedron
        coordinate, weight = gaussQuadratureTetra(GQN, FT)
        return GaussQuadratureInfoStruct{FT}(coordinate, weight, uvw)
    # 六面体
    elseif GeoS == :Hexahedron
        coordinate, weight, uvw = gaussQuadratureHexa(GQN, FT)
        return GaussQuadratureInfoStruct{FT}(coordinate, weight, uvw)
    # 其它，报错
    else
        throw("GeoS: $GeoS not supported, check the docs of this function.")            
    end
end

"""num为三角形高斯求积点数，本函数给出1, 3, 4, 7的值，默认为3, 返回值为一个坐标数组(3*num)和一个权重向量(num)"""
function gaussQuadratureTri( num::Integer = 3, FT::DataType=Float64)
    # 根据求积点数计算对应坐标权重
    if num == 1
        coordinate  =   SArray{Tuple{3, num}, FT}([1/3 1/3 1/3]) 
        weight      =   SVector{num, FT}(1.)
    elseif num == 3
        coordinate  =   SArray{Tuple{3, num}, FT}(
                            [   2/3 1/6 1/6;
                                1/6 2/3 1/6;
                                1/6 1/6 2/3])
        weight      =   SVector{num, FT}([1/3 1/3 1/3])
    elseif num == 4
        coordinate  =   SArray{Tuple{3, num}, FT}(
                            [   1/3 0.6 0.2 0.2;
                                1/3 0.2 0.6 0.2;
                                1/3 0.2 0.2 0.6])
        weight      =   SVector{num, FT}([-9/16 25/48 25/48 25/48])
    elseif num == 6
        coordinate  =   SArray{Tuple{3, num}, FT}(
                            [   0.10810301  0.44594849  0.44594849  0.81684757  0.09157621  0.09157621;
                                0.44594849  0.10810301  0.44594849  0.09157621  0.81684757  0.09157621;
                                0.44594849  0.44594849  0.10810301  0.09157621  0.09157621  0.81684757])
        weight      =   SVector{num, FT}([0.22338158 0.22338158 0.22338158 0.10995174 0.10995174 0.10995174])
        
    elseif num == 7
        a1 = 0.059715871790; b1 = 0.470142064105
        a2 = 0.797426985353; b2 = 0.101286507323
        coordinate  =   SArray{Tuple{3, num}, FT}(
                            [   1/3 a1 b1 b1 a2 b2 b2;
                                1/3 b1 a1 b1 b2 a2 b2;
                                1/3 b1 b1 a1 b2 b2 a2])
        weight      =   SVector{num, FT}([  9/40   0.132394152788506 0.132394152788506 0.132394152788506 0.125939180544827 0.125939180544827 0.125939180544827])
    else
        throw("三角形高斯求积点数 $num 目前还不支持，请输入 [1, 3, 4, 6, 7 ] 中的一个.")
        
    end # if

    return coordinate, weight
    
end # function

"""num为四面体高斯求积点数，本函数给出1, 4, 5, 11的值，默认为4, 返回值为一个坐标数组(4*num)和一个权重向量(num)"""
function gaussQuadratureTetra( num::Integer = 4, FT::DataType=Float64)

    if num == 1
        coordinate  =   SArray{Tuple{4, num}, FT}([1/4 1/4 1/4 1/4])
        weight      =   SVector{num, FT}(1.)
    elseif num == 4
        a = 0.58541020; b = 0.13819660
        coordinate  =   SArray{Tuple{4, num}, FT}(
                            [   a b b b;
                                b a b b;
                                b b a b;
                                b b b a])
        weight      =   SVector{num, FT}([1/4 1/4 1/4 1/4])
    elseif num == 5
        a = 1/2; b = 1/6
        coordinate  =   SArray{Tuple{4, num}, FT}(
                            [   1/4   a  b  b  b;
                                1/4   b  a  b  b;
                                1/4   b  b  a  b;
                                1/4   b  b  b  a])
        weight      =   SVector{num, FT}([-4/5 9/20 9/20 9/20 9/20])
        
    elseif num == 11
        a1 = 0.714285714286;b1 = 0.095238095238
        a2 = 0.399403576167;b2 = 0.100596423833
        coordinate  =   SArray{Tuple{4, num}, FT}(
                            [   1/4   a1 b1 b1 b1 a2 b2 b2 a2 b2 a2;
                                1/4   b1 a1 b1 b1 a2 a2 b2 b2 a2 b2;
                                1/4   b1 b1 a1 b1 b2 a2 a2 b2 b2 a2;
                                1/4   b1 b1 b1 a1 b2 b2 a2 a2 a2 b2])
        weight      =   SVector{num, FT}([- 0.078933333333     0.0457333333333    0.0457333333333    0.0457333333333    0.0457333333333    0.1493333333333    0.1493333333333    0.1493333333333    0.1493333333333    0.1493333333333    0.1493333333333])
    else
        throw("四面体高斯求积点数 $num 目前还不支持，请输入 [1 4 5 11] 中的一个.")
        
    end # if

    return coordinate, weight
    
end # function

"""
num 为四边形高斯求积点数，本函数给出 1, 4, 9 ... i² 的值，默认为4, 返回值为一个坐标数组(4*num)和一个权重向量(num)
"""
function gaussQuadratureQuad( num::Integer = 4, FT::DataType=Float64)

    if typeof(num) <: Integer
        # 判断输入是否有效
        n       =   num^(1/2)
        nint    =   round(Int, n)
        !isapprox(nint, n, atol=1e-8) && throw("Quadrangle Quadrature points num fails.")
        # 计算
        XGL, WGL    =   gausslegendre(nint)
        X   =   0.5 .+ 0.5 .* XGL
        W   =   0.5 .* WGL
        # uv 坐标
        uvs =   zeros(FT, 2, num)
        Ws  =   zeros(FT, num)
        for j in 1:nint, i in 1:nint
            linindex = (j-1)*nint + i
            uvs[1, linindex] = X[i]
            uvs[2, linindex] = X[j]
            Ws[linindex] = W[j] * W[i] 
        end
        # 计算求积点信息局部坐标信息（ P₍ ）
        coordinate =   zeros(4, num)
        for ii in 1:num
            # uv 坐标
            u, v  = view(uvs, :, ii)
            # 1-4
            coordinate[1, ii]   =   (1-u)*(1-v)
            coordinate[2, ii]   =      u *(1-v)
            coordinate[3, ii]   =      u *   v 
            coordinate[4, ii]   =   (1-u)*   v 
        end
        
        weight = SVector{num}(Ws)
        return  SMatrix{4, num, FT}(coordinate), weight, SMatrix{2, num, FT}(uvs)
    else
        throw("四边形高斯求积点数 $num 目前还不支持，请输入 i² 中的一值.")
    end # if

end # function

"""
num 为六面体高斯求积点数，本函数给出 1, 8, 27... i³ 的值，默认为8, 返回值为一个坐标数组(4*num)和一个权重向量(num)
"""
function gaussQuadratureHexa( num::Integer = 8, FT::DataType=Float64)

    if typeof(num) <: Integer
        # 判断输入是否有效
        n       =   num^(1/3)
        nint    =   round(Int, n)
        !isapprox(nint, n, atol=1e-8) && throw("Hexahedron Quadrature points num fails.")
        # 计算
        XGL, WGL    =   gausslegendre(nint)
        X   =   0.5 .+ 0.5 .* XGL
        W   =   0.5 .* WGL
        # uvw坐标
        uvws    =   zeros(FT, 3, num)
        Ws      =   zeros(FT, num)
        for k in 1:nint, j in 1:nint, i in 1:nint
            linindex = (k-1)*nint^2 + (j-1)*nint + i
            uvws[1, linindex] = X[i]
            uvws[2, linindex] = X[j]
            uvws[3, linindex] = X[k]
            Ws[linindex] = W[k] * W[j] * W[i] 
        end

        # 计算求积点信息局部坐标信息（ P₍ ）
        coordinate =   zeros(8, num)
        for ii in 1:num
            # uvw 坐标
            u, v, w = view(uvws, :, ii)
            # 1-8
            coordinate[1, ii]   =   (1-u)*(1-v)*(1-w)
            coordinate[2, ii]   =      u *(1-v)*(1-w)
            coordinate[3, ii]   =      u *   v *(1-w)
            coordinate[4, ii]   =   (1-u)*   v *(1-w)
            coordinate[5, ii]   =   (1-u)*(1-v)*   w
            coordinate[6, ii]   =      u *(1-v)*   w
            coordinate[7, ii]   =      u *   v *   w
            coordinate[8, ii]   =   (1-u)*   v *   w
        end
        
        weight = SVector{num, FT}(Ws)
        return SMatrix{8, num, FT}(coordinate), weight, SMatrix{3, num, FT}(uvws)
    else
        throw("六面体高斯求积点数 $num 目前还不支持，请输入 i³ 中的一值.")
    end # if

end # function


end # module GaussQuadrature4Geo