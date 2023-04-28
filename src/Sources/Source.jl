## 本文件定义源类型和一些相关的计算函数
# 源抽象类
abstract type ExcitingSource end
abstract type AntennaType <:ExcitingSource end

const ExcitingSources = Union{ExcitingSource, AbstractVector{ExcitingSource}}

# 各种类型的源
# 平面波
include("Planewave.jl")
# 磁偶极子
include("MagneticDipole.jl")
# 天线阵
include("AntettaArray.jl")

##