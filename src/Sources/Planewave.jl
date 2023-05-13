
"""
    PlaneWave{FT<:Real}<:ExcitingSource

平面波源：
```
θ   ::FT            球坐标角度θ
ϕ   ::FT            球坐标角度ϕ
α   ::FT            波极化方向相对于 θhat_source 绕K̂_v旋转的角度
f   ::FT            波频率
V   ::FT            波激励幅度
E_v ::SVec3D{FT}    入射波电场极化矢量
k̂   ::SVec3D{FT}    入射波波矢向量
```

默认构造函数：
```julia
PlaneWave{FT}(θ::FT, ϕ::FT, α::FT, V::FT = one(FT))
```
"""
struct PlaneWave{FT<:Real}<:ExcitingSource

    θ   ::FT
    ϕ   ::FT
    α   ::FT
    V   ::FT

    E_v ::SVec3D{FT}
    k̂   ::SVec3D{FT}
    
    function PlaneWave{FT}(θ::FT, ϕ::FT, α::FT, V::FT = one(FT)) where{FT<:Real}
        # 角度信息
        θϕ      =   θϕInfo{FT}(θ, ϕ)
        # 球坐标的单位向量矩阵（即直角坐标 → 球坐标转换矩阵）
        sphereUniVecs   =   cart2sphereMat(θϕ)
        # 计算电场
        E_v     =   -cos(α)*sphereUniVecs[2,:] + sin(α)*sphereUniVecs[3,:]
        # 若幅值不为1则修改
        (V != 1) && (E_v *= V)
        # 计算波矢矢量
        k̂       =   -sphereUniVecs[1,:]

        return new(θ, ϕ, α, V, E_v, k̂)
    end

end

@doc """
    PlaneWave{FT}(θ, ϕ, α, V = one(FT)) where {FT}

类型自动转换的 `PlaneWave` 构造函数。
"""
function PlaneWave{FT}(θ, ϕ, α, V = one(FT)) where {FT}
    θ::FT = θ
    ϕ::FT = ϕ
    α::FT = α
    V::FT = V
    PlaneWave{FT}(θ::FT, ϕ::FT, α::FT, V::FT)
end
PlaneWave(args...) = PlaneWave{Precision.FT}(args...) 

@doc """
    sourceEfield(plw::PlaneWave, r)

计算平面波 `plw` 在全局坐标下给定位置 `r` 处的电场。
"""
function sourceEfield(plw::PlaneWave, r)
    plw.E_v*exp(-Params.JK_0*(plw.k̂ ⋅ r))
end

@doc """
    sourceHfield(plw::PlaneWave, r)

计算平面波 `plw` 在全局坐标下给定位置 `r` 处的磁场。
"""
function sourceHfield(plw::PlaneWave, r)
    1/η_0*(plw.k̂ × sourceEfield(plw, r))
end
