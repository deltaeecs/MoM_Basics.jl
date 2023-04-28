using Base: AbstractCartesianIndex
## 导入包含的网格类型文件

abstract type VSCellType end

## 面元
abstract type SurfaceCellType{IT<:Integer, FT<:Real} <:VSCellType end
# 三角形网格类
include("Triangles.jl")

## 体元
abstract type VolumeCellType{IT<:Integer, FT<:Real, CT<:Complex{FT}} <:VSCellType end
# 四面体
include("Tetrahedras.jl")
# 六面体
include("Hexahedras.jl")


"""
此函数用于设置三角形网格介电常数。目前为空以方便多重派发。
"""
function setGeosPermittivity!(geosInfo::AbstractVector{VT}, εᵣ::CT = 1.0(1+0im)) where {VT<:TriangleInfo, CT<:Complex}

    nothing

end


"""
此函数用于设置四面体、六面体的介电常数、介质对比度，修改此函数以得到对应的数据
"""
function setGeosPermittivity!(geosInfo::AbstractVector{VT}, εᵣ::CT = 1.0(1+0im)) where {VT<:VSCellType, CT<:Complex}
    # 循环设置
    @threads for ig in eachindex(geosInfo)
        geo   =   geosInfo[ig]
        geo.ε =   εᵣ*ε_0
        geo.κ =   1 - 1/εᵣ
    end
    try
        setδκ!(geosInfo)
    catch
        @warn "设置网格间介质差异失败"
    end

    # 更新保存的参数信息
    # 更新保存的参数信息
    open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
        @printf f "%21s\n" "介质信息"
        @printf f "%-5s %28s\n" "εᵣ" εᵣ
    end

    nothing

end

"""
此函数用于设置四面体、六面体的介电常数、介质对比度，修改此函数以得到对应的数据
"""
function setGeosPermittivity!(geosInfo::AbstractVector{VT}, εᵣs::T) where {VT<:VSCellType, T<:AbstractVector}
    # 循环设置
    @threads for ig in eachindex(geosInfo)
        geo   =   geosInfo[ig]
        εᵣ    =   εᵣs[ig]
        geo.ε =   εᵣ*ε_0
        geo.κ =   1 - 1/εᵣ
    end
    try
        setδκ!(geosInfo)
    catch
        @warn "设置网格间介质差异失败"
    end

    # 更新保存的参数信息
    open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
        @printf f "%21s\n" "介质信息"
        @printf f "%-5s %28s\n" "εᵣ" mean(εᵣs)
    end

    nothing

end

"""
此函数用于设置四面体、六面体的介电常数、介质对比度，修改此函数以得到对应的数据
"""
function setGeosPermittivity!(geosInfoV::AbstractVector{VT}, εᵣ::CT = DielectricInitialε) where {VT<:AbstractVector, CT<:Complex}
    # 循环设置
    for geosInfo in geosInfoV
        setGeosPermittivity!(geosInfo, εᵣ)
    end
    nothing

end


##