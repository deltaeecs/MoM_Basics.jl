"""
阵列天线
"""
abstract type AbstractAntennaArray <:ExcitingSource end



include("TaylorWins.jl")
using .TaylorWins

"""
阵列天线
机扫阵 (mechanically scanned array, MSA)，相控阵 (Phased Array)
orient 采用的是欧拉角，参考['eulerRotationMat'](@ref)
注意阵列初始指向由提供的天线单元合成，作为阵列只提供指向旋转
"""
Base.@kwdef mutable struct AntennaArray{FT<:Real, AT, N} <: AbstractAntennaArray
    center::MVec3D{FT}  =   zero(MVec3D{FT})
    orient::MVec3D{FT}  =   zero(MVec3D{FT})
    size::NTuple{N, Int}=   Tuple(fill(0,N))
    dgap::NTuple{N, FT} =   Tuple(fill(Params.λ_0/2,N))
    antennas    ::AbstractArray{AT, N}
end

"""
初始化
"""
function AntennaArray(antennas::AbstractArray{AT, N}; center = zero(MVec3D{Precison.FT}), orient = MVec3D{Precison.FT}(0, 0, 1)) where {AT, N}
    N > 3 && throw("不支持3维以上阵列！")
    arysize =   size(antennas)
    dgap    =   Tuple([distance(eachslice(antennas; dims = i)) for i in 1:N])
    AntennaArray{Precison.FT, AT, N}(center, orient, arysize, dgap, antennas)
end

distance(arys)  =   norm(first(arys)[2].centergb - first(arys)[1].centergb)

"""
提供快捷的阵列构建函数
注意此处输入的阵列、单元指向必须为指定的欧拉角 (ZXZ) ['eulerRotationMat'](@ref)!
"""
function antennaArray(arysize, aryorient, dgap = Tuple(fill(Params.λ_0/2, length(arysize)));
    sourceConstructer, sourceT, sourceorientlc, orientunit, coefftype = :uniform, arycenter = zero(MVec3D{Precision.FT}))
    
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

    AntennaArray{FT, sourceT, length(arysize)}(arycenter, aryorient, arysize, dgap, antennas)

end

#TODO
"""
更新指向(阵列机械旋转)
ary 阵列
aryorient 阵列指向（ZXZ）
sourceorientlc 源在阵列坐标的指向欧拉角（ZXZ）
orientunit 指向欧拉角单位
"""
function update_orient!(ary::AT; aryorient, sourceorientlc, orientunit) where {AT<:AbstractAntennaArray}
    
    FT = Precision.FT
    # 将单元指向从 阵列局部坐标 转换到 全局坐标
    # 从阵列坐标到全局坐标的旋转矩阵
    aryl2gRot       =   eulerRotationMat(aryorient..., orientunit)
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

"""
更新指向相位
"""
function update_phase!(ary::AT, phasef) where {AT<:AbstractAntennaArray}
    nothing
end


"""
设置差阵列, 将一半的天线相位置反
"""
function setdiffArray!(ary, dim)
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
给定点的远场电场计算函数,\n
输入参数:\n
r_observe       ::  Vec3D{FT}, 要计算的位置(局部坐标下),\n
antennaarray    ::  Vector{ST}, 包含多个源实例的向量,\n
必选参数\n
r_coortype  ::  Char, 输入坐标类型: 直角坐标为符号： :C，球坐标： :S，默认为球坐标输入。\n
返回值：\n
电场值       ::  SVec3D{FT}\n
"""
function sourceFarEfield(ary::AT, r̂θϕ::r̂θϕInfo{FT}) where {FT<:Real, AT<:AbstractAntennaArray}
    θϕ = r̂θϕ.θϕ
    sinθ, cosθ, sinϕ, cosϕ = θϕ.sinθ, θϕ.cosθ, θϕ.sinϕ, θϕ.cosϕ
    re = zero(MVector{2, Complex{FT}})
    for source in ary.antennas
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
function sourceEfield(ary::AT, rvec::AbstractVector{FT}) where {FT<:Real,  AT<:AbstractAntennaArray}
    re = zero(MVec3D{Complex{FT}})
    for source in ary.antennas
        re .+= sourceEfield(source, rvec)
    end
    return re
end # function