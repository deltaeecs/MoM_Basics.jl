"""
阵列天线抽象类
"""
abstract type AbstractAntennaArray <:ExcitingSource end

include("TaylorWins.jl")
@reexport using .TaylorWins

"""
    AntennaArray{FT<:Real, AT, N} <: AbstractAntennaArray

机扫阵列天线 (mechanically scanned array, MSA)、相控阵 (Phased Array)阵列天线
orient 采用的是欧拉角，参考[`eulerRotationMat`](@ref)
注意阵列初始指向由提供的天线单元合成，作为阵列只提供指向旋转。
"""
Base.@kwdef mutable struct AntennaArray{FT<:Real, AT, N} <: AbstractAntennaArray
    center::MVec3D{FT}  =   zero(MVec3D{FT})
    orient::MVec3D{FT}  =   zero(MVec3D{FT})
    size::NTuple{N, Int}=   Tuple(fill(0,N))
    dgap::NTuple{N, FT} =   Tuple(fill(Params.λ_0/2,N))
    antennas::AbstractArray{AT, N}
    l2gRot::MMatrix{3, 3, FT, 9}
end

"""
    AntennaArray(antennas::AbstractArray{AT, N}; center = zero(MVec3D{Precison.FT}), orient = (0., 0., 0.)) where {AT, N}

初始化阵列天线，默认为方阵。
"""
function AntennaArray(antennas::AbstractArray{AT, N}; center = zero(MVec3D{Precison.FT}), orient = (0., 0., 0.), orientunit = :rad) where {AT, N}
    N > 3 && throw("不支持3维以上阵列！")
    arysize =   size(antennas)
    dgap    =   Tuple([distance(eachslice(antennas; dims = i)) for i in 1:N])
    l2gRot  =   eulerRotationMat(orient..., orientunit)
    AntennaArray{Precison.FT, AT, N}(center, orient, arysize, dgap, antennas, l2gRot)
end

distance(arys)  =   norm(first(arys)[2].centergb - first(arys)[1].centergb)

"""
    antennaArray(arysize, aryorient, dgap = Tuple(fill(Params.λ_0/2, length(arysize)));
    sourceConstructer, sourceT, sourceorientlc[, orientunit=:rad, coefftype = :uniform, arycenter = zero(MVec3D{Precision.FT})])

提供快捷的阵列构建函数。注意此处输入的阵列、单元指向必须为指定的欧拉角 (ZXZ) [`eulerRotationMat`](@ref)。
"""
function antennaArray(arysize, aryorient, dgap = Tuple(fill(Params.λ_0/2, length(arysize)));
    sourceConstructer, sourceT, sourceorientlc, orientunit=:rad, coefftype = :uniform, arycenter = zero(MVec3D{Precision.FT}))
    
    FT = Precision.FT
    # 将单元指向从 阵列局部坐标 转换到 全局坐标
    # 从 阵列坐标 到 全局坐标
    aryl2gRot       =   eulerRotationMat(aryorient..., orientunit)
    # 从单元局部坐标到全局坐标的旋转矩阵
    sourcel2gRot    =   aryl2gRot * eulerRotationMat(sourceorientlc..., orientunit)
    # 对应的欧拉角
    soruceorient    =   eulerRMat2αβγ(sourcel2gRot)
    # 创建偶极子阵列
    antennas     =   Array{sourceT}(undef, arysize)
    
    ## 设置具体的阵列
    # 阵列大小
    nx  =   arysize[1]; 
    ny  =   try arysize[2]; catch; 1; end
    nz  =   try arysize[3]; catch; 1; end

    dgapx::FT   =   dgap[1];
    dgapy::FT   =   try dgap[2]; catch; 0.; end
    dgapz::FT   =   try dgap[3]; catch; 0.; end
    # 阵列各个方向系数
    if coefftype == :uniform
        Imlxs   =   ones(nx)
        Imlys   =   ones(ny)
        Imlzs   =   ones(nz)
    elseif coefftype == :taylor
        Imlxs   =   nx == 1 ? ones(nx) : taylorwin(nx, 10)
        Imlys   =   ny == 1 ? ones(ny) : taylorwin(ny, 10)
        Imlzs   =   nz == 1 ? ones(nz) : taylorwin(nz, 10)
    else
        throw("目前只支持 均匀阵(:uniform) 和 泰勒阵(:taylor), 可自行设置！")
    end
    # 循环设置天线单元
    for k in 1:nz, j in 1:ny, i in 1:nx
        Iml::FT    =   Imlxs[i]*Imlys[j]*Imlzs[k]
        # 初始化
        antenna     =   sourceConstructer(;orient = MVec3D{FT}(soruceorient), Iml = Iml);
        # 单元在阵列中的位置
        centerlc    =   SVec3D{FT}((i - nx / 2 - 0.5)*dgapx, (j - ny / 2 - 0.5)*dgapy, (k - nz / 2 - 0.5)*dgapz)
        antenna.centerlc  .=  centerlc
        # 全局位置
        antenna.centergb  .=  aryl2gRot*centerlc + arycenter
        # 设置结果
        antennas[(k-1)*ny*nx + (j-1)*nx + i] = antenna
        
    end

    AntennaArray{FT, sourceT, length(arysize)}(arycenter, aryorient, arysize, dgap, antennas, aryl2gRot)

end

"""
    update_orient!(ary::AT; aryorient, sourceorientlc[, orientunit = :rad]) where {AT<:AbstractAntennaArray}

通过机械旋转更新天线阵列 `ary` 的阵列指向为 `aryorient`，天线单元指向为 `sourceorientlc`，指向角单位为 `orientunit` 。
"""
function update_orient!(ary::AT; aryorient, sourceorientlc, orientunit = :rad) where {AT<:AbstractAntennaArray}
    
    FT = Precision.FT
    if orientunit == :rad
        ary.orient .= aryorient
    elseif orientunit == :deg
        ary.orient .= aryorient .*(π/180)
    end
    # 将单元指向从 阵列局部坐标 转换到 全局坐标
    # 从阵列坐标到全局坐标的旋转矩阵
    aryl2gRot       =   eulerRotationMat(aryorient..., orientunit)
    ary.l2gRot     .=   aryl2gRot
    # 从单元局部坐标到全局坐标的旋转矩阵
    sourcel2gRot    =   aryl2gRot * eulerRotationMat(sourceorientlc..., orientunit)
    # 对应的欧拉角
    sourceorient    =   eulerRMat2αβγ(sourcel2gRot)
    # 循环设置天线单元
    for antenna in ary.antennas
        # 更新指向
        update_orient!(antenna, sourceorient, orientunit)
        # 单元在阵列中的位置
        centerlc  =  antenna.centerlc
        # 全局位置
        antenna.centergb  .=  aryl2gRot*centerlc + ary.center
    end

    nothing

end

#TODO
"""
    update_phase!(ary::AT, phase) where {AT<:AbstractAntennaArray}

更新指向相位
"""
function update_phase!(ary::AT, phase) where {AT<:AbstractAntennaArray}
    nothing
end


"""
    setdiffArray!(ary[, dim = 1])

将阵列天线 `ary` 在 `dim` 方向一半单元设置为反相位，从而实现差方向图。
"""
function setdiffArray!(ary, dim = 1)
    dim > 2 && throw("只支持一、二维阵列")
    for i in 1:ary.size[1], j in 1:ary.size[2]
        antenna =   ary.antennas[i, j]
        if dim  == 1
            (i <= ary.size[1] ÷ 2) && add_phase!(antenna, π)
        else dim == 2
            (j <= ary.size[2] ÷ 2) && add_phase!(antenna, π)
        end
    end

    nothing

end

"""
    sourceLocalFarEfield(ary::AT, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real, AT<:AbstractAntennaArray}

计算天线阵列 `ary` 在阵列坐标下给定方向 `r̂θϕ` 的远场电场。
"""
function sourceLocalFarEfield(ary::AT, r̂θϕ::r̂θϕInfo{FT}) where {FT<:AbstractFloat, AT<:AbstractAntennaArray}
    θϕ = r̂θϕ.θϕ
    sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ
    # r_xyz转换到全局坐标下
    r_xyzgb =   localrvec2Global(r̂θϕ.r̂, ary.l2gRot)
    # r̂θϕ 转换到全局坐标下
    r̂θϕgb   =   r̂θϕInfo(r_xyzgb)

    # 全局坐标下电场 θϕ 分量结果
    reθϕ = zero(MVector{2, Complex{FT}})
    for source in ary.antennas
        reθϕ .+= sourceFarEfield(source, r̂θϕgb) * exp(Params.JK_0 * (sinθ * cosϕ * source.centergb[1] + sinθ * sinϕ * source.centergb[2] + cosθ * source.centergb[3]))
    end

    # 转换到局部坐标下
    re = zero(MVec3D{Complex{FT}})
    re = (r̂θϕgb.θhat * reθϕ[1]) .+ (r̂θϕgb.ϕhat * reθϕ[2])
    temp = ary.l2gRot' * re
    # 局部坐标下的 θϕ 分量
    reθϕ[1] = r̂θϕgb.θhat ⋅ temp
    reθϕ[2] = r̂θϕgb.ϕhat ⋅ temp

    return reθϕ
end # function

"""
    sourceFarEfield(ary::AT, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real, AT<:AbstractAntennaArray}

计算天线阵列 `ary` 在全局坐标下给定方向 `r̂θϕ` 的远场电场。
"""
function sourceFarEfield(ary::AT, r̂θϕ::r̂θϕInfo{FT}) where {FT<:AbstractFloat, AT<:AbstractAntennaArray}
    θϕ = r̂θϕ.θϕ
    sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ

    # 电场结果
    reθϕ = zero(MVector{2, Complex{FT}})
    for source in ary.antennas
        reθϕ .+= sourceFarEfield(source, r̂θϕ) * exp(Params.JK_0 * (sinθ * cosϕ * source.centergb[1] + sinθ * sinϕ * source.centergb[2] + cosθ * source.centergb[3]))
    end
    return reθϕ
end # function

"""
    sourceEfield(ary::AT, rvec::AbstractVector{FT}) where {FT<:Real,  AT<:AbstractAntennaArray}

计算天线阵列 `ary` 在全局坐标下给定位置 `rvec` 处的电场。
"""
function sourceEfield(ary::AT, rvec::AbstractVector{FT}; r_coortype::Symbol=:C) where {FT<:Real,  AT<:AbstractAntennaArray}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_coortype == :C ?  (r_xyz = rvec) : (
    r_coortype == :S ?  (r_xyz = sphere2cart(rvec)) : throw("Error, 输入坐标只支持直角坐标：'C'，球坐标：'S'"))
    
    # 全局坐标电场
    re = zero(MVec3D{Complex{FT}})
    for source in ary.antennas
        re .+= sourceEfield(source, r_xyz)
    end

    # 阵列局部坐标电场
    return re
end # function

"""
    sourceLocalEfield(ary::AT, rvec::AbstractVector{FT}) where {FT<:Real,  AT<:AbstractAntennaArray}

计算天线阵列 `ary` 在阵列局部坐标下给定位置 `rvec` 处的电场。
"""
function sourceLocalEfield(ary::AT, rvec::AbstractVector{FT};  r_coortype::Symbol=:C) where {FT<:Real,  AT<:AbstractAntennaArray}
    # 笛卡尔坐标下的坐标值    r_xyz
    # 判断坐标类型并计算赋值
    r_coortype == :C ?  (r_xyz = rvec) : (
    r_coortype == :S ?  (r_xyz = sphere2cart(rvec)) : throw("Error, 输入坐标只支持直角坐标：'C'，球坐标：'S'"))
    # r_xyz转换到全局坐标下
    r_xyzgb =   localrvec2Global(r_xyz, ary.l2gRot)
    # 局部坐标电场
    re = zero(MVec3D{Complex{FT}})
    for source in ary.antennas
        re .+= sourceEfield(source, r_xyzgb)
    end
    # 转换到全局坐标
    return localrvec2Global(re, ary.l2gRot)
end # function