module MoM_Basics

using Dates, ProgressMeter, Printf
using StaticArrays, OffsetArrays, SparseArrays
using LinearAlgebra, FastGaussQuadrature, Statistics
using .Threads, ThreadsX
using Rotations
using Meshes

export  Vec3D, SVec3D, MVec3D, random_rhat,
        ∠Info, θϕInfo, r̂func, θhatfunc, ϕhatfunc, r̂θϕInfo,
        θϕInfofromCart, nodes2Poles,
        globalObs2LocalObs, localObs2GlobalObs,
        globalrvec2Local, localrvec2Global,
        dist, greenfunc,
        VecCart2SphereMat, cart2sphereMat, cart2sphere,
        eulerZunit, eulerRotationMat, eulerRMat2αβγ,
        sincmath, sphere2cart,
        IntDtype, Precision, 
        Params, modiParams!, setPrecision!,
        inputBasicParameters, saveSimulationParams,
        C_0, μ_0, ε_0, η_0, div4π, ηdiv16π,
        GaussQuadratureInfo,
        SimulationParamsType, SimulationParams, modiSimulationParams!,
        VSBFTypes, updateVSBFTypes!, updateVSBFTParams!,
        VSCellType, SurfaceCellType, VolumeCellType,
        TriangleMesh, TriangleInfo, GQPNTri, GQPNTriSglr, 
        TriGQInfo, TriGQInfoSglr, getGQPTri, getGQPTriSglr, 
        setTricoor!, setTriParam!, getGQPTri, getGQPTriSglr,
        TetrahedraInfo, GQPNTetra, GQPNTetraSglr, Tris4Tetra,
        TetraGQInfo, TetraGQInfoSglr, setTetraCoor!, setTetraParam!,
        getTetrasInfo, setδκ!, getGQPTetra, getGQPTetraSglr, 
        HexahedraMesh, Quads4Hexa, GQPNQuad1D, GQPNQuad1DSglr, GQPNQuad1DSSglr, 
        GQPNQuad, GQPNQuadSglr, GQPNQuadSSglr, QuadGQInfo, QuadGQInfoSglr, QuadGQInfoSSglr,
        getGQPQuad, getGQPQuadSglr, getGQPQuadSSglr, Quads4Hexa,
        HexahedraInfo, GQPNHexa, GQPNHexaSglr, GQPNHexaSSglr,
        HexaGQInfo, HexaGQInfoSglr, HexaGQInfoSSglr,
        setHexaCoor!, setHexaParam!, getHexasInfo,
        getGQPHexa, getGQPHexaSglr, getGQPHexaSSglr, 
        setGeosPermittivity!,
        MeshDataType, getNodeTriTetraFekoNas, getNodeTriTetraHexaNas, 
        getCellsBFs, getCellsFromFileName, getdatNodeElementParam,
        getCellsBFsFromFileName, getBFsFromMeshData, getMeshData, getConnectionMatrix,
        RWG, PWC, SWG, RBF, LinearBasisFunction, ConstBasisFunction, BasisFunctionType,
        ExcitingSource, AntennaType, ExcitingSources,
        PlaneWave, sourceEfield, sourceHfield,
        MagneticDipole, update_phase!, add_phase!, update_orient!,
        sourceFarEfield, radiationIntegralL0, radiationIntensityU_m,
        radiationPower, radiationDirectionCoeff,
        AbstractAntennaArray, taylorwin, AntennaArray, distance,
        antennaArray, setdiffArray!,
        timer, memory, @clock, show_memory_time


## 记录程序运行内存占用情况
include("Recorder.jl")

## 各部分函数
# 网格元高斯求积点、权重计算函数
include("GaussQuadrature4Geos.jl")
using .GaussQuadrature4Geo

# 一些重要的要用到的基础类定义
include("BasicStuff.jl")
include("CoorTrans.jl")
# 一些有用的函数
include("UsefulFunctions.jl")

# 参数
include("ParametersSet.jl")
# 参数输入输出
include("Inputs.jl")

# 处理网格、基函数相关
include("MeshAndBFs.jl")

## 源信息
include("Sources/Source.jl")

end
