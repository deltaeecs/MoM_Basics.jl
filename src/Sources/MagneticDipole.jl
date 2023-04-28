"""
磁偶极子天线类型\n
属性：\n
id      ::IT    编号\n
Iml     ::CT    磁流线值\n
V       ::FT    磁流线幅值
phase   ::FT    相位\n
orient  ::MVec3D{FT}  指向欧拉角\n
centerlc::MVec3D{FT}  局部坐标下的中心位置\n
centergb::MVec3D{FT}  全局坐标下的中心位置\n
l2gMat  ::MMatrix{3, 3, FT, 9}  局部坐标到全局坐标的旋转变换矩阵\n
"""
mutable struct MagneticDipole{FT<: Real}<:AntennaType
    id      ::Integer
    Iml     ::Complex{FT}
    V       ::FT
    phase   ::FT
    orient  ::MVec3D{FT}
    centerlc::MVec3D{FT}
    centergb::MVec3D{FT}
    l2gRot  ::MMatrix{3, 3, FT, 9}
end

"""
磁偶极子天线默认构造函数
id      ::IT    编号\n
Iml     ::FT    磁流线值\n
phase   ::FT    相位（输入弧度(rad)单位）\n
orient  ::MVec3D{FT}  指向欧拉角\n
centerlc::MVec3D{FT}  局部坐标下的中心位置\n
centergb::MVec3D{FT}  全局坐标下的中心位置\n
I0S     ::FT    电流环幅值 
"""
function MagneticDipole{FT}(
    id      ::Int32 =   zero(Int32);
    Iml     ::CT    =   zero(CT),
    phase   ::FT    =   zero(FT),
    orient  ::MVec3D{FT}    =   zero(MVec3D{FT}),
    centerlc::MVec3D{FT}    =   zero(MVec3D{FT}),
    centergb::MVec3D{FT}    =   zero(MVec3D{FT}),
    I0S     ::FT    =   zero(FT), unit = :rad) where{FT <: Real, CT <: Complex{FT}}

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

"""
精度可变的 MagneticDipole 构造函数
"""
function MagneticDipole(id = 0; Iml = 0., phase = 0., orient = zero(MVec3D{Float32}), centerlc = zero(MVec3D{Float32}), centergb = zero(MVec3D{Float32}))
    FT  =   Precision.FT
    CT  =   Complex{FT}
    id      ::Int32 =   id
    Iml     ::FT    =   Iml
    phase   ::FT    =   phase
    orient  ::MVec3D{FT}    =   orient
    centerlc::MVec3D{FT}    =   centerlc
    centergb::MVec3D{FT}    =   centergb
    ImlCT   =   CT(Iml*cos(phase), Iml*sin(phase))

    MagneticDipole{FT}(id; Iml = ImlCT, phase=phase, orient=orient, centerlc=centerlc, centergb=centergb)
end

"""
设置相位
"""
function update_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}
    # 当前幅值
    V   =   md.V
    # 更新磁流大小（带相位）
    md.Iml   =   CT(V*cos(phase), V*sin(phase))

    nothing
end

"""
附加相位
"""
function add_phase!(md::MagneticDipole{FT}, phase) where {FT <: Real}
    # 当前幅值
    Iml     =   md.Iml
    # 更新磁流大小（带相位）
    md.Iml  =   Iml*exp(phase*im)

    nothing
end


"""
更新指向
md:: 偶极子
orient::指向欧拉角（ZXZ）
unit::单位（度(:deg)或弧度(rad)）
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
给定点的电场计算函数,\n
输入参数:\n
r_observe   ::  Vec3D{FT}, 要计算的位置(局部坐标下),\n
magnetDipole::MagneticDipole{FT}, 要计算的磁流源实例,\n
必选参数\n
r_coortype  ::  Char, 输入坐标类型: 直角坐标为符号： :C，球坐标： :S，默认为球坐标输入。\n
返回值：\n
电场值       ::  SVec3D{FT}\n
"""
function sourceEfield(magnetDipole::MagneticDipole{FT}, r_observe::Vec3D{FT};  r_coortype::Symbol=:C) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_coortype == :C ?  (r_xyz = r_observe) : (
    r_coortype == :S ?  (r_xyz = sphere2cart(r_observe)) : throw("Error, 输入坐标只支持直角坐标：'C'，球坐标：'S'"))
    
    # r_xyz转换到 天线 局部坐标下
    r_xyzlc =   globalrvec2Local(r_xyz, magnetDipole.l2gRot, magnetDipole.centergb)
    
    # 计算 天线 局部坐标（直角）的球坐标三角函数信息
    sinθ, _, sinϕ, cosϕ  =   θϕInfofromCart(r_xyzlc)
    # 局部坐标ϕ方向向量
    ϕhat    =   SVec3D{FT}(-sinϕ, cosϕ, zero(FT))

    # 场点距离缝隙天线距离
    R_mn    =   norm(r_xyzlc)
    divR_mn =   1/R_mn
    # 局部坐标下电场的xyz三个分量
    Eveclc  =  -magnetDipole.Iml*greenfunc(R_mn)*div4π*(Params.JK_0 + divR_mn)*sinθ*ϕhat
    # 转换到全局坐标
    Evecgb  =   localrvec2Global(Eveclc, magnetDipole.l2gRot)

    return Evecgb
end # function


"""
给定点的电场计算函数,\n
输入参数:\n
r_observe   ::  Vec3D{FT}, 要计算的位置(局部坐标下),\n
magnetDipole::MagneticDipole{FT}, 要计算的磁流源实例,\n
必选参数\n
r_coortype  ::  Char, 输入坐标类型: 直角坐标为符号： :C，球坐标： :S，默认为球坐标输入。\n
返回值：\n
电场值       ::  SVec3D{FT}\n
"""
function sourceFarEfield(magnetDipole::MagneticDipole{FT}, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real}
    CT  =   Complex{FT}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_xyz   =   r̂θϕ.r̂
    # r_xyz转换到 天线 局部坐标下
    r_xyzlc =   globalrvec2Local(r_xyz, magnetDipole.l2gRot)
    
    # 计算 天线 局部坐标（直角）的球坐标三角函数信息
    sinθ, cosθ, sinϕ, cosϕ  =   θϕInfofromCart(r_xyzlc)
    # 局部坐标ϕ方向向量
    ϕhat    =   SVec3D{FT}(-sinϕ, cosϕ, zero(FT))

    # 局部坐标下电场的xyz三个分量
    Eveclc  =  -magnetDipole.Iml*div4π*(Params.JK_0)*sinθ*ϕhat
    # 转换到全局坐标
    Evecgb  =   localrvec2Global(Eveclc, magnetDipole.l2gRot)
    # 提取两个分量
    re  =   @SVector [r̂θϕ.θhat ⋅ Evecgb, r̂θϕ.ϕhat ⋅ Evecgb]

    return re
end # function


"""
给定点的远场电场计算函数,\n
输入参数:\n
r_observe   ::  Vec3D{FT}, 要计算的位置(局部坐标下),\n
sources     ::  Vector{ST}, 包含多个源实例的向量,\n
必选参数\n
r_coortype  ::  Char, 输入坐标类型: 直角坐标为符号： :C，球坐标： :S，默认为球坐标输入。\n
返回值：\n
电场值       ::  SVec3D{FT}\n
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
给定点的电场计算函数,\n
输入参数:\n
r_observe   ::  Vec3D{FT}, 要计算的位置(局部坐标下),\n
sources     ::  Vector{ST}, 包含多个源实例的向量,\n
必选参数\n
r_coortype  ::  Char, 输入坐标类型: 直角坐标为符号： :C，球坐标： :S，默认为球坐标输入。\n
返回值：\n
电场值       ::  SVec3D{FT}\n
"""
function sourceEfield(sources::Vector{ST}, rvec::AbstractVector{FT}) where {FT<:Real, ST<:ExcitingSource}
    re = zero(MVec3D{Complex{FT}})
    for source in sources
        re .+= sourceEfield(source, rvec)
    end
    return re
end # function


"""
磁偶极子的远场辐射积分计算函数，注意θϕ为偶极子的局部坐标
L(θ, ϕ)
输入参数:
θϕ      ::θϕInfo{FT},   要计算的空间角度信息，对远场采用球坐标计算
magnetDipole::MagneticDipole{FT}, 要计算的磁流源实例,\n
返回值：
电场值  ::MVector{2, FT}，在θ、ϕ方向的计算结果
"""
function radiationIntegralL0(θϕ::θϕInfo{FT}, magnetDipole::MagneticDipole{FT}) where {FT<:Real}
    # 结果
    L_θ =  -magnetDipole.Iml*θϕ.sinθ
    L_ϕ =   zero(FT)
    # 结果数组
    return SVector{2, FT}(L_θ, L_ϕ)

end # function


"""
磁流源的辐射强度计算函数
U_m(θ, ϕ) = Y_0/(8λ_0²)*(|L_θ|² + |L_ϕ|²)
"""
function radiationIntensityU_m(θϕ::θϕInfo{FT}, magnetDipole::MagneticDipole{FT}) where {FT<:Real}
    # 先计算磁流源的辐射积分
    local L = radiationIntegralL0(θϕ, magnetDipole)
    # 再计算辐射强度
    return Y_0/(8*Params.λ_0^2)*(abs2(L[1]) + abs2(L[2]))
end # function

"""
辐射功率计算函数
P_rad = ∫∫  U(θ, ϕ)sinθ  dθdϕ
对磁偶极子可直接在源缝表面积分：
P_rad = ∫∫ |E(r)|²/(2η₀) dxdy
    """
radiationPower(magnetDipole::MagneticDipole{FT}) where {FT<:Real} = Y_0*π/(3*Params.λ_0^2)*abs2(magnetDipole.Iml)

"""
方向性系数
D_m(θ, ϕ) = 4π U_m(θ, ϕ)/P_rad
"""
function radiationDirectionCoeff(θϕ::θϕInfo{FT}, magnetDipole::MagneticDipole{FT}) where {FT<:Real}
    # 先计算磁流源的辐射强度
    local U =   radiationIntensityU_m(θϕ, magnetDipole)
    
    # P_rad 辐射功率计算
    P_rad   =   radiationPower(magnetDipole)

    # 计算辐射方向性系数并返回
    return 4π*U/P_rad
end # function
