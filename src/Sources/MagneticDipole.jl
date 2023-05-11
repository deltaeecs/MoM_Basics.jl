"""
    MagneticDipole{FT<: Real}<:AntennaType

磁偶极子天线类型。
```
id      ::Integer               编号
Iml     ::Complex{FT}           磁流线值
V       ::FT                    磁流线幅值
phase   ::FT                    相位
orient  ::MVec3D{FT}            指向欧拉角
centerlc::MVec3D{FT}            局部坐标下的中心位置
centergb::MVec3D{FT}            全局坐标下的中心位置
l2gRot  ::MMatrix{3, 3, FT, 9}  局部坐标到全局坐标的旋转变换矩阵
```
"""
mutable struct MagneticDipole{FT<: Real}<:AntennaType
    id      ::Int
    Iml     ::Complex{FT}
    V       ::FT
    phase   ::FT
    orient  ::MVec3D{FT}
    centerlc::MVec3D{FT}
    centergb::MVec3D{FT}
    l2gRot  ::MMatrix{3, 3, FT, 9}
end

"""
```
    MagneticDipole{FT}(
    id      ::Int32         =   zero(Int32);        # 编号
    Iml     ::CT            =   zero(CT),           # 磁流线值
    phase   ::FT            =   zero(FT),           # 相位（输入弧度(rad)单位）
    orient  ::MVec3D{FT}    =   zero(MVec3D{FT}),   # 指向欧拉角
    centerlc::MVec3D{FT}    =   zero(MVec3D{FT}),   # 局部坐标下的中心位置
    centergb::MVec3D{FT}    =   zero(MVec3D{FT}),   # 全局坐标下的中心位置
    I0S     ::FT            =   zero(FT),           # 电流环幅值
    unit                    =   :rad    ) 
    where{FT <: Real, CT <: Complex{FT}}
```
"""
function MagneticDipole{FT}(
    id      ::IT    =   zero(Int32);
    Iml     ::CT    =   zero(Complex{FT}),
    phase   ::FT    =   zero(FT),
    orient  ::MVec3D{FT}    =   zero(MVec3D{FT}),
    centerlc::MVec3D{FT}    =   zero(MVec3D{FT}),
    centergb::MVec3D{FT}    =   zero(MVec3D{FT}),
    I0S     ::FT    =   zero(FT), unit = :rad) where{IT<:Integer, FT <: Real, CT <: Complex{FT}}

    ImlFTtemp   =   zero(FT)
    ## 根据输入的源同意计算Iml
    if iszero(Iml) & iszero(I0S)
        # 默认为幅度为1的Iml
        ImlFTtemp = one(CT)
    elseif iszero(Iml) & ~iszero(I0S)
        # 输入为I0S
        ImlFTtemp = Params.ω_0*μ_0*I0S
    elseif ~iszero(Iml) & iszero(I0S)
        # 输入为Iml
        ImlFTtemp = convert(CT, Iml)
    elseif ~iszero(Iml) & ~iszero(I0S)
        # 输入有问题
        throw("只能输入一种类型的源")
    else
        # 输入有问题
        throw("检查源幅度输入，只支持浮点数。")
    end
    # 计算Iml在相位phase下的复数形式
    ImlCT   =   CT(ImlFTtemp*cos(phase), ImlFTtemp*sin(phase))
    # 计算实数幅值 V
    V       =   abs(ImlCT)
    # 计算局部坐标到全局坐标的旋转变换矩阵
    l2gRot  =   eulerRotationMat(orient..., unit)
    
    MagneticDipole{FT}(id, ImlCT, V, phase, orient, centerlc, centergb, l2gRot)

end
MagneticDipole(args...) = MagneticDipole{Precision.FT}(args...)

"""
    update_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}

设置磁偶极 `md` 的相位为 `phase`。
"""
function update_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}
    # 当前幅值
    V   =   md.V
    # 更新磁流大小（带相位）
    md.Iml   =   CT(V*cos(phase), V*sin(phase))

    nothing
end

"""
    add_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}

为磁偶极 `md` 附加相位 `phase`。
"""
function add_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}
    # 当前幅值
    Iml     =   md.Iml
    # 更新磁流大小（带相位）
    md.Iml  =   Iml*exp(phase*im)

    nothing
end


"""
    update_orient!(md::MagneticDipole{FT}, orient, unit = :rad) where {FT <: Real}

更新磁偶极 `md` 指向为 `orient`。
"""
function update_orient!(md::MagneticDipole{FT}, orient, unit = :rad) where {FT <: Real}
    # 更新指向
    md.orient  .=   if unit === :rad
        orient
    elseif unit === :deg
        orient ./ 180 * π
    else
        throw("角度单位设置错误")
    end
    # 更新全局坐标到本地坐标的旋转矩阵
    md.l2gRot  .=   eulerRotationMat(orient..., unit)

    nothing
end

"""
    sourceLocalEfield(md::MagneticDipole{FT}, r_observe::Vec3D{FT};  r_coortype::Symbol=:C) where {FT<:Real}

计算磁偶极 `md` 在磁偶极局部坐标给定位置 `rvec` 处的电场。
"""
function sourceLocalEfield(md::MagneticDipole{FT}, r_observe::Vec3D{FT};  r_coortype::Symbol=:C) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_coortype == :C ?  (r_xyz = r_observe) : (
    r_coortype == :S ?  (r_xyz = sphere2cart(r_observe)) : throw("Error, 输入坐标只支持直角坐标：'C'，球坐标：'S'"))
    
    # 计算 天线 局部坐标（直角）的球坐标三角函数信息
    sinθ, _, sinϕ, cosϕ  =   θϕInfofromCart(r_xyz)
    # 局部坐标ϕ方向向量
    ϕhat    =   SVec3D{FT}(-sinϕ, cosϕ, zero(FT))

    # 场点距离缝隙天线距离
    R_mn    =   norm(r_xyz)
    divR_mn =   1/R_mn
    # 局部坐标下电场的xyz三个分量
    Eveclc  =  -md.Iml*greenfunc(R_mn)*div4π*(Params.JK_0 + divR_mn)*sinθ*ϕhat

    return Eveclc
end # function

"""
    sourceEfield(md::MagneticDipole{FT}, r_observe::Vec3D{FT};  r_coortype::Symbol=:C) where {FT<:Real}

计算磁偶极 `md` 在全局坐标下给定位置 `rvec` 处的电场。
"""
function sourceEfield(md::MagneticDipole{FT}, r_observe::Vec3D{FT};  r_coortype::Symbol=:C) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_coortype == :C ?  (r_xyz = r_observe) : (
    r_coortype == :S ?  (r_xyz = sphere2cart(r_observe)) : throw("Error, 输入坐标只支持直角坐标：'C'，球坐标：'S'"))
    
    # r_xyz转换到 天线 局部坐标下
    r_xyzlc =   globalrvec2Local(r_xyz, md.l2gRot, md.centergb)
    # 局部坐标下电场的xyz三个分量
    Eveclc  =  sourceLocalEfield(md, r_xyzlc)
    # 转换到全局坐标
    Evecgb  =   localrvec2Global(Eveclc, md.l2gRot)

    return Evecgb
end # function

"""
    sourceLocalFarEfield(md::MagneticDipole{FT}, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real}

计算磁偶极 `md` 在磁偶极局部坐标下给定方向 `r̂θϕ` 的远场电场。
"""
function sourceLocalFarEfield(md::MagneticDipole, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值
    r_xyzlc =   r̂θϕ.r̂
    # θϕ
    θϕ = r̂θϕ.θϕ

    # 局部坐标ϕ方向向量
    ϕhat    =   r̂θϕ.ϕhat

    # 局部坐标下电场的 ϕ 分量
    re  =   @SVector [0, -md.Iml*div4π*(Params.JK_0)*θϕ.sinθ]

    return re
end # function


"""
    sourceFarEfield(md::MagneticDipole{FT}, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real}

计算磁偶极 `md` 在全局坐标下给定方向 `r̂θϕ` 的远场电场。
"""
function sourceFarEfield(md::MagneticDipole, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_xyz   =   r̂θϕ.r̂
    # r_xyz转换到 天线 局部坐标下
    r_xyzlc =   globalrvec2Local(r_xyz, md.l2gRot)
    
    # 计算 天线 局部坐标（直角）的球坐标三角函数信息
    sinθ, _, sinϕ, cosϕ  =   θϕInfofromCart(r_xyzlc)
    # 局部坐标ϕ方向向量
    ϕhat    =   SVec3D{FT}(-sinϕ, cosϕ, zero(FT))

    # 局部坐标下电场的xyz三个分量
    Eveclc  =  -md.Iml*div4π*(Params.JK_0)*sinθ*ϕhat
    # 转换到全局坐标
    Evecgb  =   localrvec2Global(Eveclc, md.l2gRot)
    # 提取两个分量
    re  =   @SVector [r̂θϕ.θhat ⋅ Evecgb, r̂θϕ.ϕhat ⋅ Evecgb]

    return re
end # function


"""
    sourceFarEfield(sources::Vector{ST}, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real, ST<:ExcitingSource}

计算源向量 `sources` 在全局坐标下给定方向 `r̂θϕ` 的远场电场。
"""
function sourceFarEfield(sources::Vector{ST}, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real, ST<:ExcitingSource}
    θϕ = r̂θϕ.θϕ
    sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ
    re = zero(MVector{2, Complex{FT}})
    for source in sources
        re .+= sourceFarEfield(source, r̂θϕ) * exp(Params.JK_0 * (sinθ * cosϕ * source.centergb[1] + sinθ * sinϕ * source.centergb[2] + cosθ * source.centergb[3]))
    end
    return re
end # function

"""
    sourceEfield(sources::Vector{ST}, rvec::AbstractVector{FT}) where {FT<:Real, ST<:ExcitingSource}

计算源向量 `sources` 在全局坐标下给定位置 `rvec` 处的远场电场。
"""
function sourceEfield(sources::Vector{ST}, rvec::AbstractVector{FT}) where {FT<:Real, ST<:ExcitingSource}
    re = zero(MVec3D{Complex{FT}})
    for source in sources
        re .+= sourceEfield(source, rvec)
    end
    return re
end # function


"""
    radiationIntegralL0(md::MagneticDipole, θϕ::θϕInfo{FT}) where {FT<:Real}

磁偶极子的远场辐射积分计算函数，注意 `θϕ` 为偶极子的局部坐标。
"""
function radiationIntegralL0(md::MagneticDipole, θϕ::θϕInfo{FT}) where {FT<:Real}
    # 结果
    L_θ =  -md.Iml*θϕ.sinθ
    L_ϕ =   zero(FT)
    # 结果数组
    return SVector{2, FT}(L_θ, L_ϕ)

end # function


@doc raw"""
    radiationIntensityU_m(md::MagneticDipole{FT}, θϕ::θϕInfo{FT}) where {FT<:Real}

计算磁流源的辐射强度函数
``U_m(θ, ϕ) = \frac{Y_0}{8λ_0²}(|L_θ|² + |L_ϕ|²)``。
"""
function radiationIntensityU_m(md::MagneticDipole{FT}, θϕ::θϕInfo{FT}) where {FT<:Real}
    # 先计算磁流源的辐射积分
    local L = radiationIntegralL0(md, θϕ)
    # 再计算辐射强度
    return Y_0/(8*Params.λ_0^2)*(abs2(L[1]) + abs2(L[2]))
end # function

@doc raw"""
    radiationPower(md::MagneticDipole{FT}) where {FT<:Real}

计算辐射功率。
``P_{rad} = ∫∫  U(θ, ϕ)sinθ  dθdϕ``
对磁偶极子可直接在源缝表面积分：
``P_{rad} = ∫∫ |E(r)|²/(2η₀) dxdy``
"""
radiationPower(md::MagneticDipole) = Y_0*π/(3*Params.λ_0^2)*abs2(md.Iml)

@doc raw"""
radiationDirectionCoeff(md::MagneticDipole{FT}, θϕ::θϕInfo{FT}) where {FT<:Real}

计算方向性系数：``D_m(θ, ϕ) = 4π U_m(θ, ϕ)/P_{rad}``。
"""
function radiationDirectionCoeff(md::MagneticDipole{FT}, θϕ::θϕInfo{FT}) where {FT<:Real}
    # 先计算磁流源的辐射强度
    local U =   radiationIntensityU_m(md, θϕ)
    
    # P_rad 辐射功率计算
    P_rad   =   radiationPower(md)

    # 计算辐射方向性系数并返回
    return 4π*U/P_rad
end # function
