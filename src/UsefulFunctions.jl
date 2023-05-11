@doc """
    dist(pa::AbstractVector{FT}, pb::AbstractVector{FT})::FT where {FT<:AbstractFloat}
    dist(pa::Vec3D{FT}, pb::Vec3D{FT})::FT where {FT<:AbstractFloat}
    dist(pa::Vec3D{FT})::FT where {FT<:AbstractFloat}

计算两点之间距离，比使用norm函数更高效。
"""
function dist(pa::AbstractVector{FT}, pb::AbstractVector{FT})::FT where {FT<:AbstractFloat}
    @inbounds sqrt((pa[1]-pb[1])^2 + (pa[2]-pb[2])^2 + (pa[3]-pb[3])^2)
end #function
function dist(pa::Vec3D{FT}, pb::Vec3D{FT})::FT where {FT<:AbstractFloat}
    @inbounds sqrt((pa[1]-pb[1])^2 + (pa[2]-pb[2])^2 + (pa[3]-pb[3])^2)
end #function
function dist(pa::Vec3D{FT})::FT where {FT<:AbstractFloat}
    @inbounds sqrt(abs2(pa[1]) + abs2(pa[2]) + abs2(pa[3]))
end #function

@doc """
    greenfunc(R::T) where {T<:AbstractFloat}
    greenfunc(pa::Vec3D{T}, pb::Vec3D{T}) where {T<:AbstractFloat}
    greenfunc(pa::AbstractVector{T}, pb::AbstractVector{T}) where {T<:AbstractFloat}
    greenfunc(R::T, k::T) where {T<:AbstractFloat}
    greenfunc(pa::Vec3D{T}, pb::Vec3D{T}, k::T) where {T<:AbstractFloat}

计算归一化自由空间格林函数
``g(R) =  exp^{-1im*K_0*R}/R``
"""
function greenfunc(R::T) where {T<:AbstractFloat}
    exp(-Params.JK_0*R)/R
end #function
function greenfunc(pa::Vec3D{T}, pb::Vec3D{T}) where {T<:AbstractFloat}
    # 两点距离
    local R =   dist(pa, pb)
    # 格林函数
    return greenfunc(R)
end
function greenfunc(pa::AbstractVector{T}, pb::AbstractVector{T}) where {T<:AbstractFloat}
    # 两点距离
    local R =   dist(pa, pb)
    # 格林函数
    return greenfunc(R)
end
function greenfunc(R::T, k::T) where {T<:AbstractFloat}
    exp(-im*k*R)/R
end #function
function greenfunc(pa::Vec3D{T}, pb::Vec3D{T}, k::T) where {T<:AbstractFloat}
    # 两点距离
    local R =   dist(pa, pb)
    # 格林函数
    return greenfunc(R, k)
end

"""
    acrossbcrossc(a, b, c)

更高效地计算矢量连续叉乘：``a × b × c = (c⋅a)b - (c⋅b)a``
"""
acrossbcrossc(a, b, c) = b*(c⋅a) - a*(c⋅b)

"""
    VecCart2SphereMat{T} <: FieldVector{3, Vec3D{T}}

矢量坐标系（直角 → 球）转换矩阵
```
r_hat   ::  Vec3DCart{T}    r 方向的单位向量
θ_hat   ::  Vec3DCart{T}    θ 方向的单位向量
ϕ_hat   ::  Vec3DCart{T}    ϕ 方向的单位向量
```
"""
struct VecCart2SphereMat{T} <: FieldVector{3, Vec3D{T}}
    r_hat   ::Vec3D{T} 
    θ_hat   ::Vec3D{T}
    ϕ_hat   ::Vec3D{T}
    """
    VecCart2SphereMat 默认构造函数
    """
    function VecCart2SphereMat{FT}() where{FT<:AbstractFloat}
        new(zero(VecCart2SphereMat{FT}))
    end #function VecCart2Cart
end # struct



"""
    VecCart2SphereMat{FT}(θϕ::θϕInfo{FT}) where{FT<:AbstractFloat}

`VecCart2SphereMat`在给定方向 `θϕ` 的构造函数。
"""
function VecCart2SphereMat{FT}(θϕ::θϕInfo{FT}) where{FT<:AbstractFloat}
    # 先分解角度相关信息
    local sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ

    # 结果预分配
    local re  = VecCart2SphereMat{FT}()
    # 结果计算
    re.r_hat    .=  Vec3DCart{FT}([sinθ*cosϕ  sinθ*sinϕ   cosθ])
    re.θ_hat    .=  Vec3DCart{FT}([cosθ*cosϕ  cosθ*sinϕ  -sinθ])
    re.ϕ_hat    .=  Vec3DCart{FT}([    -sinϕ       cosϕ   zero(FT)])
    return re
end # function

"""
    Base.:*(c2smat::VecCart2SphereMat{FT}, vec3D::Vec3D{T}) where {T<:Number, FT<:AbstractFloat}

重载 * 函数以计算向量在球坐标下的各个分量
"""
function Base.:*(c2smat::VecCart2SphereMat{FT}, vec3D::Vec3D{T}) where {T<:Number, FT<:AbstractFloat}
    MVec3D{T}(c2smat.r_hat ⋅ vec3D, c2smat.θ_hat ⋅ vec3D, c2smat.ϕ_hat ⋅ vec3D)
end

@doc """
    cart2sphereMat(θϕ::θϕInfo{FT}) where{FT<:Real}
    cart2sphereMat(θ::FT, ϕ::FT) where{FT<:Real}

计算矢量坐标系（直角 → 球）转换矩阵。
"""
function cart2sphereMat(θϕ::θϕInfo{FT}) where{FT<:Real}
    # 分解角度相关信息
    sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ
    # 计算转换矩阵
    mat = MMatrix{3, 3, FT}([   sinθ*cosϕ  sinθ*sinϕ   cosθ     ; 
                                cosθ*cosϕ  cosθ*sinϕ  -sinθ     ;
                                    -sinϕ       cosϕ   zero(FT)  ])
    return mat
end
function cart2sphereMat(θ::FT, ϕ::FT) where{FT<:Real}
    # 计算θϕInfo
    θϕ  =   θϕInfo{FT}(θ, ϕ)
    # 利用多重派发的另一个函数计算
    return cart2sphereMat(θϕ)
end

"""
    eulerZunit(α::FT, β::FT, γ::FT, unit::Symbol) where{FT<:Real}

计算给定欧拉角局部坐标 `z` 轴在全局坐标中的位置单位向量。
"""
function eulerZunit(α::FT, β::FT, γ::FT, unit::Symbol) where{FT<:Real}
    # 角度信息转换到弧度 rad
    α_rad   =   zero(FT)
    β_rad   =   zero(FT)
    γ_rad   =   zero(FT)
    # 计算旋转角度的三角函数信息
    if unit == :rad 
        α_rad   =   α
        β_rad   =   β
        γ_rad   =   γ
    elseif unit == :deg
        α_rad   =   α/180*π
        β_rad   =   β/180*π
        γ_rad   =   γ/180*π
    else
        throw("请确定坐标类型，取 (:rad 或 :deg)")
    end
    # 计算
    SVec3D{FT}(sin(β_rad)*sin(α_rad), -sin(β)*cos(α_rad), cos(β_rad))
end


"""
    eulerRotationMat(α::FT, β::FT, γ::FT, unit::Symbol) where{FT<:Real}

根据坐标旋转的欧拉角计算旋转矩阵, 定义旋转顺序为：
“滚动” → “俯仰” → “偏航”，
即按绕 “z轴” → “x轴” → “z轴”的顺序，分别旋转
α, β, γ 度
[Wikipedia-Euler_angles](https://en.wikipedia.org/wiki/Euler_angles)
输入：
α, β, γ, 旋转角度信息
unit: 输入角度值单位，默认为 :rad，可选 :deg
输出：
rotMat :: SMatrix{3, 3, FT}, 坐标旋转矩阵
rotMat * vec 将 vec 从局部坐标转换回全局坐标
"""
function eulerRotationMat(α::FT, β::FT, γ::FT, unit::Symbol) where{FT<:Real}
    # 角度信息转换到弧度 rad
    α_rad   =   zero(FT)
    β_rad   =   zero(FT)
    γ_rad   =   zero(FT)
    # 计算旋转角度的三角函数信息
    if unit == :rad 
        α_rad   =   α
        β_rad   =   β
        γ_rad   =   γ
    elseif unit == :deg
        α_rad   =   α/180*π
        β_rad   =   β/180*π
        γ_rad   =   γ/180*π
    else
        throw("请确定坐标类型")
    end

    # 计算并返回
    return MMatrix{3, 3}(RotZXZ{FT}(α_rad, β_rad, γ_rad))
end


"""
    eulerRotationMat(axis::Vec3D{FT}, θ::FT, unit::Symbol) where{FT<:Real}

计算天线阵按给定任意轴 `axis`, 旋转 `θ`` 角度的旋转矩阵。
"""
function eulerRotationMat(axis::Vec3D{FT}, θ::FT, unit::Symbol) where{FT<:Real}

    # 计算旋转角度的三角函数信息
    θ_rad   = if unit == :rad 
        θ
    elseif unit == :deg
        θ/180*π
    else
        throw("请确定坐标类型")
    end
    # 利用Rotations包计算并返回
    return MMatrix{3, 3}(AngleAxis{FT}(θ_rad, axis...))
end

"""
    eulerRotationMat(θ::FT, ϕ::FT, unit::Symbol) where{FT<:AbstractFloat}

计算转到给定指向 `θ, ϕ` 处的旋转矩阵，旋转一步到位，不发生自旋。
"""
function eulerRotationMat(θ::FT, ϕ::FT, unit::Symbol) where{FT<:AbstractFloat}

    θ_rad   =   zero(FT)
    ϕ_rad   =   zero(FT)
    # 计算旋转角度的三角函数信息
    if unit == :rad 
        θ_rad   =   θ
        ϕ_rad   =   ϕ
    elseif unit == :deg
        θ_rad   =   θ/180*π
        ϕ_rad   =   ϕ/180*π
    else
        throw("请确定坐标类型")
    end

    # 计算旋转轴
    axis    =  [cos(ϕ_rad + π/2), sin(ϕ_rad + π/2), 0.]

    # 利用Rotations包计算并返回
    return MMatrix{3, 3}(AngleAxis{FT}(θ_rad, axis...))
end


"""
    eulerRMat2αβγ(l2gMat)

通过局部坐标到全局坐标的旋转矩阵 `l2gMat` 计算欧拉角。
旋转矩阵按绕 “z轴” → “x轴” → “z轴”的顺序，分别旋转
`α, β, γ` 度得到，本函数计算对应的三个角度 `α, β, γ`。
"""
function eulerRMat2αβγ(l2gMat)
    # 计算局部 ẑ 在全局坐标的位置
    ẑgb =   l2gMat * [0, 0, 1.]
    # 计算 α
    α   =   atan(ẑgb[1], -ẑgb[2])
    # 计算 β
    β   =   acos(ẑgb[3])
    # 计算 x̂, ŷ 旋转  “z轴” → “x轴”  和 x̂ 完成全部旋转后到达的位置矢量
    l2gZX   =   RotZX{eltype(l2gMat)}(α, β)
    x̂zxgb   =   l2gZX*[1., 0, 0]
    ŷzxgb   =   l2gZX*[0, 1., 0]
    x̂gb     =   l2gMat * [1., 0, 0]

    # x̂gb 在 x̂zxgb, ŷzxgb 的投影
    cosγ    =   x̂gb ⋅ x̂zxgb
    sinγ    =   x̂gb ⋅ ŷzxgb
    # 计算 
    γ   =   atan(sinγ, cosγ)

    return α, β, γ
end

const divπ = 1/π
"""
    sincmath(x::T) where{T<:Number}

Julia 自带 [sinc](@ref Base.sinc) 函数计算的是归一化辛格函数:

``sinc(x)     =   sin(πx)/(πx)``

此处借用 `sinc`，定义数学领域的非归一化 sinc 函数，即计算:

``sin(x)/x``

"""
sincmath(x::T) where{T<:Number} = sinc(x*divπ)

@doc """
    sphere2cart(coor_sphere::Vec3D{T}) where T<:Real
    sphere2cart(coor_sphere...)
    sphere2cart(r::T, θϕ::θϕInfo{T}) where T<:Real

将球坐标 `coor_sphere` 转换到直角坐标。
"""
function sphere2cart(coor_sphere::Vec3D{T}) where T<:Real
    # 球坐标下的一些量
    local r         =   coor_sphere[1]
    local sinθ,cosθ =   sincos(coor_sphere[2])
    local sinϕ,cosϕ =   sincos(coor_sphere[3])
    # 计算结果并返回
    return  MVec3D{T}(r*sinθ*cosϕ, r*sinθ*sinϕ, r*cosθ)
end
function sphere2cart(coor_sphere...)
    sphere2cart(convert(SVec3D{mapreduce(typeof, promote_type, coor_sphere)}, coor_sphere))
end
function sphere2cart(r::T, θϕ::θϕInfo{T}) where T<:Real
    # 先分解角度相关信息
    local sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ
    # 计算结果并返回
    return  MSVec3D{T}(r*sinθ*cosϕ, r*sinθ*sinϕ, r*cosθ)
end

@doc """
    cart2sphere(xyz::Vararg{T, 3}) where {T}

将直角坐标 `xyz` 转换到球坐标。
"""
function cart2sphere(xyz::Vararg{T, 3}) where {T}
    # 直角坐标下的一些量
    x, y, z =   xyz
    # r
    r       =   sqrt(x^2 + y^2 + z^2)
    @assert !isapprox(r, 0; atol = eps(typeof(r)))
    # θ
    θ       =   acos(z/r)
    # ϕ
    ϕ       =   atan(y, x)
    ϕ < 0 && begin ϕ += 2π; end
    # 计算结果并返回
    return  MVec3D{typeof(r)}(r, θ, ϕ)
end
function cart2sphere(xyz::AbstractVector)
    cart2sphere(xyz...)
end

"""
    random_rhat(; FT = Precision.FT)

随机生成单位向量。
"""
function random_rhat(; FT = Precision.FT)
    x, y, z = let x::FT, y::FT, z::FT, r::FT
        while true
            x = 2rand()-1
            y = 2rand()-1
            z = 2rand()-1
            r = sqrt(x^2+y^2+z^2)
            r <=1 && break
        end
        x/r, y/r, z/r
    end
    SVec3D{FT}(x, y, z)
end