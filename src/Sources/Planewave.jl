
"""
平面波源，定义在局部坐标下，与球坐标下的单位向量存在如下关系
θhat_source =   - θhat
K̂_v         =   - r̂
包含如下信息、
θ   ::FT    ，入射球坐标角度θ
ϕ   ::FT    ，入射球坐标角度ϕ
α   ::FT    ，入射波极化方向相对于 θhat_source 绕K̂_v旋转的角度
f   ::FT    ，入射波频率
V   ::FT    ，入射波激励幅度
E_v ::SVec3D{FT}    ，入射波电场极化矢量
k̂   ::SVec3D{FT}    ，入射波波矢向量
"""
struct PlaneWave{FT<:Real}<:ExcitingSource

    θ   ::FT
    ϕ   ::FT
    α   ::FT
    V   ::FT

    E_v ::SVec3D{FT}
    k̂   ::SVec3D{FT}
    
    """
    默认构造函数
    输入
    θ   ::FT
    ϕ   ::FT
    α   ::FT
    f   ::FT
    V   ::FT
    计算
    E_v ::SVec3D{FT}
    k̂   ::SVec3D{FT}
    """
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

"""
类型自动转换的构造函数
"""
function PlaneWave{FT}(θ, ϕ, α, V = one(FT)) where {FT}
    θ::FT = θ
    ϕ::FT = ϕ
    α::FT = α
    V::FT = V
    PlaneWave{FT}(θ::FT, ϕ::FT, α::FT, V::FT)
end

"""
类型自动转换的构造函数
"""
function PlaneWave(θ, ϕ, α, V = 1)
    FT = Precision.FT
    θ::FT = θ
    ϕ::FT = ϕ
    α::FT = α
    V::FT = V
    PlaneWave{FT}(θ::FT, ϕ::FT, α::FT, V::FT)
end

"""
平面波源的空间r点的电场 E 计算函数
输入信息：
source  ::PlaneWave{FT}
r       ::Vec3D{FT}
返回值
E       ::Vec3D{Complex{FT}}
"""
@inline function sourceEfield(source::PlaneWave{FT}, r::Vec3D{FT}) where{FT<:Real}
    source.E_v*exp(-Params.JK_0*(source.k̂ ⋅ r))
end

"""
平面波源的空间r点的磁场 H 计算函数
输入信息：
source  ::PlaneWave{FT}
r       ::Vec3D{FT}
返回值
E       ::Vec3D{Complex{FT}}
"""
@inline function sourceHfield(source::PlaneWave{FT}, r::Vec3D{FT}) where{FT<:Real}
    1/η_0*(source.k̂ × sourceEfield(source, r))
end
