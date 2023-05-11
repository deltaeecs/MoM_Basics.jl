"""
    globalObs2LocalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::StaticMatrix{3,3, FT}) where {FT}

根据全局观测空间角度信息 `r̂θϕs_obs` 计算给定局部至全局坐标旋转矩阵 `l2gRot` 下局部坐标的观测空间角度信息。
"""
function globalObs2LocalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::StaticMatrix{3,3, FT}) where {FT}
    # 旋转后的局部指向向量 r̂
    r̂ObsLocal       =   [l2gRot'*r̂θϕ_obs.r̂ for r̂θϕ_obs in r̂θϕs_obs]
    # 根据r̂重建观测空间角度信息矩阵
    r̂θϕsObsLocal    =   [r̂θϕInfo{FT}(r̂) for r̂ in r̂ObsLocal]
    
    return r̂θϕsObsLocal
end

"""
    localObs2GlobalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::StaticMatrix{3,3, FT}) where {FT}

根据局部观测空间角度信息 `r̂θϕs_obs` 计算给定局部至全局坐标旋转矩阵 `l2gRot` 下全局坐标的观测空间角度信息。
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
    globalrvec2Local(rvecsglobal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}

计算全局向量组成的矩阵 `rvecsglobal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的局部坐标。
"""
function globalrvec2Local(rvecsglobal::Matrix{T}, l2gRot::StaticMatrix{3, 3, FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecslocal      =   similar(rvecsglobal)
    @threads for ri in axes(rvecsglobal, 2)
        @inbounds rvecslocal[:, ri]  =  l2gRot' * rvecs[:, ri]
    end
    return rvecslocal
end

"""
    globalrvec2Local(rvecsglobal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}

计算全局向量组成的矩阵 `rvecsglobal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的局部坐标，局部坐标的原点在全局坐标的 `r0InGlobal` 处。
"""
function globalrvec2Local(rvecsglobal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    # 预分配内存
    rvecslocal      =   similar(rvecsglobal)

    @threads for ri in axes(rvecsglobal, 2)
        @inbounds rvecslocal[:, ri]  =  l2gRot' * (rvecs[:, ri] .- r0InGlobal)
    end

    return rvecslocal
end


"""
    localrvec2Global(rvecslocal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}

计算局部向量组成的矩阵 `rvecslocal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的全局坐标。
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
    localrvec2Global(rvecslocal::Matrix{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}

计算局部向量组成的矩阵 `rvecslocal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的全局坐标，局部坐标的原点在全局坐标的 `r0InGlobal` 处。
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
    globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}

计算全局向量 `rvecglobal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的局部坐标。
"""
function globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot' * rvecglobal
end

"""
    globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}

计算全局向量 `rvecglobal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的局部坐标，局部坐标的原点在全局坐标的 `r0InGlobal` 处。
"""
function globalrvec2Local(rvecglobal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot' * (rvecglobal .- r0InGlobal)
end


"""
    localrvec2Global(rvecslocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}

计算局部向量 `rveclocal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的全局坐标。
"""
function localrvec2Global(rveclocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot * rveclocal
end

"""
    localrvec2Global(rvecslocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}

计算局部向量 `rveclocal` 在给定局部至全局坐标旋转矩阵 `l2gRot` 下的全局坐标，局部坐标的原点在全局坐标的 `r0InGlobal` 处。
"""
function localrvec2Global(rvecslocal::Vec3D{T}, l2gRot::StaticMatrix{3,3, FT}, r0InGlobal::Vec3D{FT}) where {T<:Number, FT<:Real}
    @inbounds l2gRot * rvecslocal .+ r0InGlobal
end
