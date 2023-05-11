## 本文件用于导入各种基函数相关文件

"""
基函数抽象类型
"""
abstract type BasisFunctionType end

"""
常数基抽象类
"""
abstract type ConstBasisFunction <: BasisFunctionType end

"""
线性基抽象类
"""
abstract type LinearBasisFunction <: BasisFunctionType end
# 面基函数
include("RWG.jl")

# 体基函数
include("PWC.jl")
include("SWG.jl")
include("RBF.jl")

##
