## 这里提供一些各处都用得到的变量

# 三维静态、动态向量
Vec3D{T}    =   StaticVector{3, T} where {T<:Number}
SVec3D{T}   =   SVector{3, T} where T<:Number
MVec3D{T}   =   MVector{3, T} where T<:Number

SVec3D{T}(x::Number) where T<:Number = fill(x, SVec3D{T})
MVec3D{T}(x::Number) where T<:Number = fill(x, MVec3D{T})
Vec3D{T}(x::Number)  where T<:Number = fill(x, Vec3D{T})

"""
    ∠Info{FT<:Real}

∠ 空间角度信息类型，保存以避免大量重复计算。
"""
struct ∠Info{FT<:Real}
    ∠   ::FT
    sin∠::FT
    cos∠::FT
end # struct

function ∠Info{FT}(∠::FT = zero(FT)) where {FT<:Real}
    local sin∠, cos∠ =   sincos(∠)
    return ∠Info{FT}(∠, sin∠, cos∠)
end # function

"""
    θϕInfo{FT<:Real}

θ ϕ 空间角度信息类型，保存以避免大量重复计算。
"""
struct θϕInfo{FT<:Real}
    θ   ::FT
    ϕ   ::FT
    sinθ::FT
    cosθ::FT
    sinϕ::FT
    cosϕ::FT
end # struct

@doc """
    θϕInfo{FT}(θ::FT = zero(FT), ϕ::FT = zero(FT)) where {FT<:Real}
    θϕInfo(θ::FT, ϕ::FT) where {FT<:Real}
    θϕInfo{FT}(θ::∠Info{FT}, ϕ::∠Info{FT}) where {FT<:Real}
    θϕInfo{FT}(θ::FT, ϕ::∠Info{FT}) where {FT<:Real}
    θϕInfo{FT}(θ::∠Info{FT}, ϕ::FT) where {FT<:Real}

输入角度 `θ` 和 `ϕ` 构造 `θϕInfo` 实例。
"""
function θϕInfo{FT}(θ::FT = zero(FT), ϕ::FT = zero(FT)) where {FT<:Real}
    local sinθ, cosθ =   sincos(θ)
    local sinϕ, cosϕ =   sincos(ϕ)
    θϕInfo{FT}(θ, ϕ, sinθ, cosθ, sinϕ, cosϕ)
end # function
θϕInfo(θ::FT, ϕ::FT) where {FT<:Real} = θϕInfo{FT}(θ, ϕ)
function θϕInfo{FT}(θ::∠Info{FT}, ϕ::∠Info{FT}) where {FT<:Real}
    θϕInfo{FT}(θ.∠, ϕ.∠, θ.sin∠, θ.cos∠, ϕ.sin∠, ϕ.cos∠)
end # function
function θϕInfo{FT}(θ::FT, ϕ::∠Info{FT}) where {FT<:Real}
    local sinθ, cosθ =   sincos(θ)
    θϕInfo{FT}(θ, ϕ.∠, sinθ, cosθ, ϕ.sin∠, ϕ.cos∠)
end # function
function θϕInfo{FT}(θ::∠Info{FT}, ϕ::FT) where {FT<:Real}
    local sinϕ, cosϕ =   sincos(ϕ)
    θϕInfo{FT}(θ.∠, ϕ, θ.sin∠, θ.cos∠, sinϕ, cosϕ)
end # function

"""
    θϕInfo{FT}(rvec::AbstractVector{FT}) where {FT<:Real}

输入直角坐标角度 `rvec` 构造 `θϕInfo` 实例。
"""
function θϕInfo{FT}(rvec::AbstractVector{FT}) where {FT<:Real}
    # xy平面的投影
    ρ   =   sqrt(rvec[1]^2 + rvec[2]^2)
    if ρ == 0
        sinϕ    =   zero(FT)
        cosϕ    =   rvec[3]/abs(rvec[3])
    else
        sinϕ    =   rvec[2]/ρ
        cosϕ    =   rvec[1]/ρ
    end
    # 向量长度
    r       =   norm(rvec)
    sinθ    =   ρ/r
    cosθ    =   rvec[3]/r

    θ   =  asin(ρ/r)
    ϕ   =  atan(rvec[2], rvec[1])

    return θϕInfo{FT}(θ, ϕ, sinθ, cosθ, sinϕ, cosϕ)
end # function

# 用于给定θϕ计算r̂的函数
r̂func(θϕ::θϕInfo{FT}) where{FT<:Real}       =   MVec3D{FT}([θϕ.sinθ*θϕ.cosϕ   θϕ.sinθ*θϕ.sinϕ   θϕ.cosθ])
r̂func(θ::FT, ϕ::FT) where{FT<:Real}         =   MVec3D{FT}([sin(θ)*cos(ϕ)  sin(θ)*sin(ϕ)   cos(θ)])
θhatfunc(θϕ::θϕInfo{FT}) where{FT<:Real}    =   MVec3D{FT}([θϕ.cosθ*θϕ.cosϕ   θϕ.cosθ*θϕ.sinϕ  -θϕ.sinθ])
θhatfunc(θ::FT, ϕ::FT) where{FT<:Real}      =   MVec3D{FT}([cos(θ)*cos(ϕ)  cos(θ)*sin(ϕ)  -sin(θ)])
ϕhatfunc(θϕ::θϕInfo{FT}) where{FT<:Real}    =   MVec3D{FT}([-θϕ.sinϕ    θϕ.cosϕ   zero(FT)])
ϕhatfunc(ϕ::FT) where{FT<:Real}             =   MVec3D{FT}([-sin(ϕ)     cos(ϕ)    zero(FT)])


"""
    r̂θϕInfo{FT<:Real}

r̂ θ ϕ 空间角度信息类型，保存以避免大量重复计算
"""
struct r̂θϕInfo{FT<:Real}
    r̂   ::MVec3D{FT}
    θhat::MVec3D{FT}
    ϕhat::MVec3D{FT}
    θϕ  ::θϕInfo{FT}
end # struct

"""
    r̂θϕInfo(θ::FT = zero(FT), ϕ::FT = zero(FT)) where {FT<:Real}
    r̂θϕInfo(θ::∠Info{FT}, ϕ::∠Info{FT}) where {FT<:Real}
    r̂θϕInfo(θ::FT, ϕ::∠Info{FT}) where {FT<:Real}
    r̂θϕInfo(θ::∠Info{FT}, ϕ::FT) where {FT<:Real}

r̂θϕInfo 的构造函数
输入角度 `θ` 和 `ϕ` 构造 `r̂θϕInfo` 实例。
"""
function r̂θϕInfo(θ::FT = zero(FT), ϕ::FT = zero(FT)) where {FT<:Real}
    local sinθ, cosθ =   sincos(θ)
    local sinϕ, cosϕ =   sincos(ϕ)
    θϕ  =   θϕInfo{FT}(θ, ϕ, sinθ, cosθ, sinϕ, cosϕ)
    r̂   =   r̂func(θϕ)
    θhat=   θhatfunc(θϕ)
    ϕhat=   ϕhatfunc(θϕ)
    r̂θϕInfo{FT}(r̂, θhat, ϕhat, θϕ)
end # function
function r̂θϕInfo(θ::∠Info{FT}, ϕ::∠Info{FT}) where {FT<:Real}
    θϕ  =   θϕInfo{FT}(θ.∠, ϕ.∠, θ.sin∠, θ.cos∠, ϕ.sin∠, ϕ.cos∠)
    r̂   =   r̂func(θϕ)
    θhat=   θhatfunc(θϕ)
    ϕhat=   ϕhatfunc(θϕ)
    r̂θϕInfo{FT}(r̂, θhat, ϕhat, θϕ)
end # function
function r̂θϕInfo(θ::FT, ϕ::∠Info{FT}) where {FT<:Real}
    local sinθ, cosθ =   sincos(θ)
    θϕ  =   θϕInfo{FT}(θ, ϕ.∠, sinθ, cosθ, ϕ.sin∠, ϕ.cos∠)
    r̂   =   r̂func(θϕ)
    θhat=   θhatfunc(θϕ)
    ϕhat=   ϕhatfunc(θϕ)
    r̂θϕInfo{FT}(r̂, θhat, ϕhat, θϕ)
end # function
function r̂θϕInfo(θ::∠Info{FT}, ϕ::FT) where {FT<:Real}
    local sinϕ, cosϕ =   sincos(ϕ)
    θϕ  =   θϕInfo{FT}(θ.∠, ϕ, θ.sin∠, θ.cos∠, sinϕ, coϕ)
    r̂   =   r̂func(θϕ)
    θhat=   θhatfunc(θϕ)
    ϕhat=   ϕhatfunc(θϕ)
    r̂θϕInfo{FT}(r̂, θhat, ϕhat, θϕ)
end # function

"""
    r̂θϕInfo{FT}(rvec::AbstractVector{FT}) where {FT<:Real}

输入直角坐标角度 `rvec` 构造 `r̂θϕInfo` 实例。
"""
function r̂θϕInfo(rvec::AbstractVector{FT}) where {FT<:Real}
    r̂   =   rvec ./ norm(rvec)
    θϕ      =   θϕInfo{FT}(r̂)
    θhat    =   θhatfunc(θϕ)
    ϕhat    =   ϕhatfunc(θϕ)
    r̂θϕInfo{FT}(r̂, θhat, ϕhat, θϕ)   
end # function

"""
    θϕInfofromCart(rvec::Vec3D{FT}) where {FT<:Real}

从直角坐标 `rvec` 计算三角函数：``sinθ, cosθ, sinϕ, cosϕ``
"""
function θϕInfofromCart(rvec::AbstractVector{FT}) where {FT<:Real}
    # xy平面的投影
    ρ   =   sqrt(rvec[1]^2 + rvec[2]^2)
    if ρ == 0
        sinϕ    =   zero(FT)
        cosϕ    =   one(FT)
    else
        sinϕ    =   rvec[2]/ρ
        cosϕ    =   rvec[1]/ρ
    end
    # 向量长度
    r       =   norm(rvec)
    sinθ    =   ρ/r
    cosθ    =   rvec[3]/r

    return sinθ, cosθ, sinϕ, cosϕ
end # function


"""
    nodes2Poles(nodes::Matrix{FT}) where {FT}

将球面散点 `nodes` 转化为角度信息实例 `r̂θϕInfo` 数组。
"""
function nodes2Poles(nodes::Matrix{FT}) where {FT}
    # 点数
    nr  =   size(nodes, 2)
    # 初始化
    r̂sθsϕs  =   Vector{r̂θϕInfo{FT}}(undef, nr)

    # 循环计算
    for ii in 1:nr
        r̂sθsϕs[ii] = r̂θϕInfo{FT}(view(nodes, :, ii))
    end

    return r̂sθsϕs

end


##