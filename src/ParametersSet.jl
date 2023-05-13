# 设置算法为单精度以节省内存
const IntDtype      =   Int64
# const FloatDtype    =   Float32
using Base: @kwdef
# const ComplexDtype  =   Complex{FloatDtype}

"""
创建类型用于控制程序精度
"""
mutable struct PrecisionT
    FT::DataType
    CT::DataType
end

PrecisionT(FT::T) where {T<:DataType} =   PrecisionT(FT, Complex{FT})

"""
创建实例控制精度
"""
Precision   =   PrecisionT(Float64)

"""
    setPrecision!(FT::Type{T}) where {T<:Union{Float32, Float64}}

设置仿真精度为 `FT`。
"""
function setPrecision!(FT::Type{T}) where {T<:Union{Float32, Float64}}
    Precision.FT = FT
    Precision.CT = Complex{FT}
    nothing
end

# const IntUnion      =   Union{Int32, Int64}
# const FloatUnion    =   Union{Float32, Float64}
# const ComplexUnion  =   Union{ComplexF32, ComplexF64}


# 真空光速
const C_0           =   299792458.0
# 真空磁导率
const μ_0           =   4e-7π
# 真空介电常数
const ε_0           =   1/(μ_0*C_0^2)
# 真空波阻抗
const η_0           =   sqrt(μ_0/ε_0)
# 真空导纳
const Y_0           =   1/η_0
# 1/4π
const div4π         =   1/4π
# η_0/(16π)
const ηdiv16π       =   η_0/(16π)

"""
ParamsType{FT<:AbstractFloat, CT<:Complex{FT}}

创建参数类型以方便保存仿真参数并控制精度：
```
frequency   ::FT        频率
ω_0         ::FT        角频率
K_0         ::FT        波数
λ_0         ::CT        波长
Jη_0divK    ::CT        常数
JK_0        ::CT        常数
k²          ::FT        常数
divJω       ::CT        常数
mk²div9     ::FT        常数
mJη_0div4πK ::CT        常数
C4divk²     ::FT        常数
JKηdiv16π   ::CT        常数
div9Jω      ::CT        常数
Rsglr       ::FT        处理奇异性的相对距离阈值
resultDir   ::String    存放结果的位置
```
"""
mutable struct ParamsType{FT<:AbstractFloat, CT<:Complex{FT}}
    frequency   ::FT
    ω_0         ::FT
    K_0         ::FT
    λ_0         ::FT
    JKη_0       ::CT
    Jη_0divK    ::CT
    JK_0        ::CT
    k²          ::FT
    divJω       ::CT
    mk²div9     ::FT
    mJη_0div4πK ::CT
    C4divk²     ::FT
    JKηdiv16π   ::CT
    div9Jω      ::CT
    Rsglr       ::FT
    CFIEα       ::FT
end

@doc """
    ParamsType{FT, CT}(frequency::FT) where{FT<:AbstractFloat, CT<:Complex{FT}}
    ParamsType(frequency::FT) where {FT<:AbstractFloat}
    ParamsType{FT}(frequency) where {FT<:AbstractFloat}

输入频率创建仿真参数实例。
"""
function ParamsType{FT, CT}(frequency::FT) where{FT<:AbstractFloat, CT<:Complex{FT}}
    # 角频率
    ω_0         =   2π*frequency
    # 真空波数
    K_0         =   ω_0/C_0
    # 波长
    λ_0         =   C_0/frequency
    # 一些用到的常量
    JKη_0       =   im*K_0*η_0
    Jη_0divK    =   im*η_0/K_0
    JK_0        =   im*K_0
    k²          =   K_0^2
    divJω       =   1/(im*ω_0)
    mk²div9     =   -k²/9
    mJη_0div4πK =   -Jη_0divK/4π
    C4divk²     =   4/k²
    JKηdiv16π   =   1im*K_0*η_0/16π
    div9Jω      =   divJω/9
    # 处理奇异性的相对距离阈值
    Rsglr       =   0.15*λ_0
    CFIEα       =   0.6

    return ParamsType{FT, Complex{FT}}( frequency, ω_0, K_0, λ_0, JKη_0, Jη_0divK, JK_0,
                                        k², divJω, mk²div9, mJη_0div4πK, C4divk², 
                                        JKηdiv16π, div9Jω, Rsglr, CFIEα)

end
ParamsType(frequency::FT) where {FT<:AbstractFloat} =  ParamsType{FT, Complex{FT}}(frequency) 
ParamsType{FT}(frequency) where {FT<:AbstractFloat} =  ParamsType{FT, Complex{FT}}(convert(FT, frequency)) 

"""
保存仿真参数的实例。
"""
const Params  =   ParamsType{Precision.FT}(5e8)

"""
    modiParams!(   frequency, ω_0, K_0, λ_0, JKη_0, Jη_0divK, JK_0,
                        k², divJω, mk²div9, mJη_0div4πK, 
                        C4divk², JKηdiv16π, div9Jω, Rsglr, CFIEα)

输入所有参数修改仿真参数的实例。
"""
function modiParams!(   frequency, ω_0, K_0, λ_0, JKη_0, Jη_0divK, JK_0,
                        k², divJω, mk²div9, mJη_0div4πK, 
                        C4divk², JKηdiv16π, div9Jω, Rsglr, CFIEα)
    Params.frequency    =   frequency  
    Params.ω_0          =   ω_0        
    Params.K_0          =   K_0        
    Params.λ_0          =   λ_0
    Params.JKη_0        =   JKη_0     
    Params.Jη_0divK     =   Jη_0divK   
    Params.JK_0         =   JK_0       
    Params.k²           =   k²         
    Params.divJω        =   divJω     
    Params.mk²div9      =   mk²div9    
    Params.mJη_0div4πK  =   mJη_0div4πK
    Params.C4divk²      =   C4divk²  
    Params.JKηdiv16π    =   JKηdiv16π  
    Params.div9Jω       =   div9Jω
    Params.Rsglr        =   Rsglr
    Params.CFIEα        =   CFIEα

    # 初始化时间和内存计数
    initialize_time_and_mem()

    BLAS.set_num_threads(Threads.nthreads())

    nothing
end

"""
    SimulationParamsType

非数值仿真参数信息：
```
resultDir   ::String        结果文件夹路径
ieT         ::Symbol        积分方程类型，包括 EFIE, MFIE, CFIE等
meshfilename::String        网格文件名称
meshunit    ::Symbol        网格文件单位
SHOWIMAGE   ::Bool          根是否要在前端显示图片
discreteVar ::String        离散的体未知量类型，支持位移电流 `"D"` 或等效电流 `"J"`
sbfT        ::Symbol        面基函数类型，目前仅支持 `:RWG`
vbfT        ::Symbol        体基函数类型，目前支持 `:SWG, :RBF, :PWC`
```
"""
mutable struct SimulationParamsType
    resultDir   ::String
    ieT         ::Symbol
    meshfilename::String
    meshunit    ::Symbol
    SHOWIMAGE   ::Bool
    discreteVar ::String
    sbfT        ::Symbol
    vbfT        ::Symbol
end

"""
    SimulationParamsType(;ieT::Symbol=:EFIE, meshfilename::String = "", meshunit::Symbol = :mm, 
    SHOWIMAGE = true, discreteVar = "D", sbfT = :RWG, vbfT = :nothing)

`SimulationParamsType` 的默认构造函数。
"""
function SimulationParamsType(;ieT::Symbol=:EFIE, meshfilename::String = "", meshunit::Symbol = :mm, 
    SHOWIMAGE = true, discreteVar = "D", sbfT = :RWG, vbfT = :nothing)
    # 保存结果的数组
    resultDir   =   "results/"*"$(Date(now()))/$(hour(now())).$(minute(now())) $(Params.frequency/1e9)GHz/"
    SimulationParamsType(resultDir, ieT, meshfilename, meshunit, SHOWIMAGE, discreteVar, sbfT, vbfT)
end

"""
非数值仿真参数实例
"""
SimulationParams    =   SimulationParamsType()

"""
    modiSimulationParams!(;ieT::Symbol=SimulationParams.ieT, 
    meshfilename::String = SimulationParams.meshfilename, 
    meshunit = SimulationParams.meshunit,
    SHOWIMAGE = SimulationParams.SHOWIMAGE,
    discreteVar = SimulationParams.discreteVar
    )

ieT         ::Symbol, 积分方程类型，包括 EFIE, MFIE, CFIE等
"""
function modiSimulationParams!(;ieT::Symbol=SimulationParams.ieT, 
    meshfilename::String = SimulationParams.meshfilename, 
    meshunit    = SimulationParams.meshunit,
    SHOWIMAGE   = SimulationParams.SHOWIMAGE,
    discreteVar = SimulationParams.discreteVar,
    sbfT        = SimulationParams.sbfT,
    vbfT        = SimulationParams.vbfT
    )
    # 保存结果的数组
    resultDir   =   "results/"*"$(Date(now()))/$(hour(now())).$(minute(now())) $(Params.frequency/1e9)GHz/"
    SimulationParams.resultDir  =   resultDir
    SimulationParams.ieT        =   ieT
    SimulationParams.meshfilename   =   meshfilename
    SimulationParams.meshunit   =   meshunit
    SimulationParams.SHOWIMAGE  =   SHOWIMAGE
    SimulationParams.discreteVar=   discreteVar
    SimulationParams.sbfT       =   sbfT
    SimulationParams.vbfT       =   vbfT

    saveSimulationParams()
    nothing
end


"""
1-3 的循环
"""
const Vec3IdxCircle =   SVector{5}([1, 2, 3, 1, 2])

"""
三角形在构建时构成第i个边的两个点为三角形中的除了第i个点的两个点，以下为索引构成第i个边的第一个点（计算边向量被 - 掉）时采用的对应点
"""
const EDGEVmINTriVsID   =   SVec3D{Int}(2, 3, 1)
"""
三角形在构建时构成第i个边的两个点为三角形中的除了第i个点的两个点，以下为索引构成第i个边的第二个点（计算边向量时用于减去第二个点 掉）时采用的对应点
"""
const EDGEVpINTriVsID   =   SVec3D{Int}(3, 1, 2)
