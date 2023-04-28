"""
根据全局观测空间角度信息计算给定坐标旋转矩阵下局部坐标的观测空间角度信息
输入
r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, 全局坐标观测信息
l2gRot::StaticMatrix{3,3, FT}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
"""
function globalObs2LocalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::StaticMatrix{3,3, FT}) where {FT}
    # 旋转后的局部指向向量 r̂
    r̂ObsLocal       =   [l2gRot'*r̂θϕ_obs.r̂ for r̂θϕ_obs in r̂θϕs_obs]
    # 根据r̂重建观测空间角度信息矩阵
    r̂θϕsObsLocal    =   [r̂θϕInfo{FT}(r̂) for r̂ in r̂ObsLocal]
    
    return r̂θϕsObsLocal
end

"""
根据局部观测空间角度信息计算给定坐标旋转矩阵下全局坐标的观测空间角度信息
输入
r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, 局部坐标观测信息
l2gRot::StaticMatrix{3,3, FT}, 局部全局坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
"""
function localObs2GlobalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::StaticMatrix{3,3, FT}) where {FT}
    # 旋转后的局部指向向量 r̂, θhat, ϕhat
    r̂ObsGlobal      =   [l2gRot*r̂θϕ_obs.r̂    for r̂θϕ_obs in r̂θϕs_obs]
    θhatObsGlobal   =   [l2gRot*r̂θϕ_obs.θhat for r̂θϕ_obs in r̂θϕs_obs]
    ϕhatObsGlobal   =   [l2gRot*r̂θϕ_obs.ϕhat for r̂θϕ_obs in r̂θϕs_obs]

    # 根据r̂重建观测空间角度信息矩阵
    r̂θϕsObsGlobal   =   [r̂θϕInfo{FT}(r̂ObsGlobal[ii], θhatObsGlobal[ii], ϕhatObsGlobal[ii]) for ii in eachindex(r̂ObsGlobal)]
    r̂θϕsObsGlobal   =   reshape(r̂θϕsObsGlobal, size(r̂θϕs_obs))
    return r̂θϕsObsGlobal
end



"""
计算全局向量在给定坐标旋转矩阵下的局部坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
"""
function globalrvec2Local(rvecglobal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecslocal      =   similar(rvecglobal)
    @threads for ri in 1:size(rvecglobal, 2)
        @inbounds rvecslocal[:, ri]  =  l2gRot' * rvecs[:, ri]
    end
    return rvecslocal
end

"""
计算全局点在给定坐标旋转矩阵下的局部坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
r0InGlobal::Vec3D{FT}局部坐标原点在全局坐标中的位置
"""
function globalrvec2Local(rvecglobal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecslocal      =   similar(rvecglobal)

    @threads for ri in 1:size(rvecglobal, 2)
        @inbounds rvecslocal[:, ri]  =  l2gRot' * (rvecs[:, ri] .- r0InGlobal)
    end

    return rvecslocal
end


"""
计算局部向量在给定坐标旋转矩阵下的全局坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 局部坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵
"""
function localrvec2Global(rvecslocal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecsglobal      =   similar(rvecslocal)
    @threads for ri in 1:size(rvecslocal, 2)
        @inbounds rvecsglobal[:, ri]  =  l2gRot * rvecs[:, ri]
    end
    return rvecsglobal
end

"""
计算全局点在给定坐标旋转矩阵下的局部坐标
输入
rvecslocal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
r0InGlobal::Vec3D{FT}局部坐标原点在全局坐标中的位置
"""
function localrvec2Global(rvecslocal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecsglobal      =   similar(rvecslocal)

    @threads for ri in 1:size(rvecslocal, 2)
        @inbounds rvecsglobal[:, ri]  =  l2gRot * rvecs[:, ri] .+ r0InGlobal
    end

    return rvecsglobal
end

"""
计算全局向量在给定坐标旋转矩阵下的局部坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
"""
@inline function globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot' * rvecglobal
end

"""
计算全局点在给定坐标旋转矩阵下的局部坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
r0InGlobal::Vec3D{FT}局部坐标原点在全局坐标中的位置
"""
@inline function globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot' * (rvecglobal .- r0InGlobal)
end


"""
计算局部向量在给定坐标旋转矩阵下的全局坐标
输入
rvecglobal::Matrix{r̂θϕInfo{FT}}, 局部坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵
"""
@inline function localrvec2Global(rvecslocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot * rvecslocal
end

"""
计算全局点在给定坐标旋转矩阵下的局部坐标
输入
rvecslocal::Matrix{r̂θϕInfo{FT}}, 全局坐标、
l2gRot ::MMatrix{3,3,F}, 局部坐标到全局坐标的旋转矩阵，因此计算时用其转置计算全局到局部
r0InGlobal::Vec3D{FT}局部坐标原点在全局坐标中的位置
"""
@inline function localrvec2Global(rvecslocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot * rvecslocal .+ r0InGlobal
end
