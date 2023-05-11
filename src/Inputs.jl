"""
    inputBasicParameters(;frequency::FT = 1e8, ieT::Symbol = :EFIE, CFIEα::FT = 0.6,
    meshfilename::String = SimulationParams.meshfilename) where {FT<:AbstractFloat}

输入频率参数 `frequency`，修改其它仿真参数的函数；
积分方程类型参数 `ieT`，修改计算过程中采用的积分方程；
CFIE混合系数 `CFIEα`、网格文件名 `meshfilename`。
"""
function inputBasicParameters(;frequency::FT = Params.frequency, ieT::Symbol = :EFIE, CFIEα::FT = 0.6,
    meshfilename::String = SimulationParams.meshfilename) where {FT<:AbstractFloat}
    # 设置 BLAS 的线程防止冲突
    @info "CEM_MoMs start with $(nthreads()) threads."
    BLAS.set_num_threads(nthreads())
    frequency::Precision.FT   =   frequency
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
    Rsglr       =   0.15λ_0
    CFIEα       =   CFIEα

    # 更改相关参数
    modiParams!(frequency, ω_0, K_0, λ_0, JKη_0, Jη_0divK, JK_0,
                k², divJω, mk²div9, mJη_0div4πK, 
                C4divk², JKηdiv16π, div9Jω, Rsglr, CFIEα)
    modiSimulationParams!(ieT = ieT, meshfilename = meshfilename)

    return nothing
end

"""
    saveSimulationParams(;meshfilename::String = SimulationParams.meshfilename, 
    sbfT::Symbol = SimulationParams.sbfT, vbfT::Symbol = SimulationParams.vbfT)

保存仿真参数到结果文件中。
"""
function saveSimulationParams(;meshfilename::String = SimulationParams.meshfilename, 
    sbfT::Symbol = SimulationParams.sbfT, vbfT::Symbol = SimulationParams.vbfT,
    targetfile = "InputArgs.txt")
    # 更新积分方程类型
    SimulationParams.sbfT = sbfT
    SimulationParams.vbfT = vbfT
    # 创建结果目录
    ~isdir(SimulationParams.resultDir) && mkpath(SimulationParams.resultDir)

    # 将输入信息保存到结果目录
    open(joinpath(SimulationParams.resultDir, targetfile), "a+")  do f
        @printf f "%19s\n" "仿真参数"
        @printf f "%-9s %24s\n" "输入文件" meshfilename
        @printf f "%-9s %24s\n" "仿真精度" Precision.FT
        @printf f "%-9s %24i\n" "线程数量" nthreads()
        @printf f "%-13s %16.4f GHz\n" "频率" Params.frequency/1e9
        @printf f "%-19s %14s\n" "积分方程类型" SimulationParams.ieT
    end

    return nothing

end