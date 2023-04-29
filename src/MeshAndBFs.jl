# 根据网格中三角形\四面体的数量，建立RWG、PWC、SWG、RBF基函数
# 网格文件处理函数
# include("MeshFileReader.jl")
# using .MeshFileReader
include("MeshProcess/MeshProcess.jl")
include("BasicVSCellType/CellTypes.jl")
include("BasisFunctions/BFs.jl")

mutable struct VSBFTstruct
    sbfType::DataType
    vbfType::DataType
end

"""
创建字典保存本次仿真时的面、体基函数类型
"""
const VSBFTypes =   VSBFTstruct(BasisFunctionType, BasisFunctionType)

"""
更新体、面基函数类型字典
"""
function updateVSBFTypes!(;sbfType = BasisFunctionType, vbfType = BasisFunctionType)
    VSBFTypes.sbfType  =   sbfType
    VSBFTypes.vbfType  =   vbfType
    nothing
end

function updateVSBFTParams!(;sbfT = :nothing, vbfT = :nothing)
    # 精度
    FloatDtype  =   Precision.FT

    # 更新面、体基函数类型
    sbfType     = begin
        sbfT == :nothing ? BasisFunctionType : begin
                sbfT == :RWG ? RWG{IntDtype, FloatDtype} : throw("面网格与基函数类型不匹配！")
            end
    end

    vbfType     = begin
        vbfT == :nothing ? BasisFunctionType : begin
            if vbfT == :SWG
                SWG{IntDtype, FloatDtype}
            elseif vbfT == :RBF
                RBF{IntDtype, FloatDtype}
            elseif vbfT == :PWC
                PWC{IntDtype, FloatDtype}
            else
                throw("体网格与基函数类型不匹配！")
            end
        end
    end

    updateVSBFTypes!(;sbfType = sbfType, vbfType = vbfType)

    nothing
end

"""
根据网格元（如三角形、四边形、四面体、六面体等）获取基函数信息用于快速计算当前单元采用的基函数
"""
getBFTfromCellT(::Type{T}) where {T<:TriangleInfo}      =   VSBFTypes.sbfType
getBFTfromCellT(::Type{T}) where {T<:VolumeCellType}    =   VSBFTypes.vbfType

"""
根据读取的网格数据生成网格元 + 基函数信息
返回值 
"""
function getCellsBFs(meshData, vbfT)

    ngeo    =   0
    nbf     =   0

    nbfs    =   Int[]
    bfTs    =   Symbol[]
    # 有三角形则创建三角形、基函数信息
    meshData.trinum     > 0     && begin
        # 计算三角形信息
        print("Constructing RWG basis function...\t")
        trianglemeshData    =   TriangleMesh(meshData.trinum, meshData.node, meshData.triangles)
        trisInfo, stbfsInfo   =   getTriangleInfo(trianglemeshData)
        # 三角形数量
        ntri    =   length(trisInfo)
        # rwg基函数数量
        nsbf    =   length(stbfsInfo)
        # 简单检查网格是否封闭，否则必须使用EFIE
        SimulationParams.ieT == :CFIE && begin
            # 单连通封闭面网格 应满足 基函数数量是三角形的 3/2 倍
            3ntri != 2nsbf && begin 
                SimulationParams.ieT = :EFIE
                @warn "网格不封闭，不可采用 CFIE, 因此恢复使用 EFIE"
            end
        end

        ngeo   +=   ntri
        nbf    +=   nsbf

        push!(nbfs, nsbf)
        push!(bfTs, :RWG)
        println("Done!")
        println("共得到 $ntri 个三角形， $nsbf 个 RWG 基函数。")
    end #begin

    # 有四面体则创建四面体、基函数信息
    meshData.tetranum   > 0     && begin
        # 开始构造四面体信息、基函数信息
        print("Constructing VIE basis function based on tetrahedras...\t")
        tetrameshData   =   TetrahedraMesh( meshData.tetranum,   meshData.node,  meshData.tetrahedras)
        if (vbfT === :PWC) || (vbfT === :SWG)
            tetrasInfo, vtbfsInfo =   getTetrasInfo(tetrameshData, vbfT)
        else
            throw("体积分基函数类型出错或与网格类型不匹配，请检查！")
        end
        # 设置介质体介电参数
        setGeosPermittivity!(tetrasInfo)
        # 三角形数量
        ntetra   =   length(tetrasInfo)
        # 基函数数量
        nvtbf    =   length(vtbfsInfo)
        # 混合网格时更新 id 信息
        @threads for geo in tetrasInfo
            geo.tetraID  +=  ngeo
            geo.inBfsID .+=  nbf
        end
        # 混合网格时更新基函数信息
        @threads for bf in vtbfsInfo
            bf.bfID      +=  nbf
            for ii in 1:length(bf.inGeo)
                bf.inGeo[ii] == 0 && continue
                if vbfT === :SWG
                    bf.inGeo[ii] +=  ngeo
                elseif vbfT === :PWC
                    bf.inGeo    +=  ngeo
                end
            end
        end
        # 几何信息更新
        ngeo    +=   ntetra
        nbf     +=   nvtbf
        push!(nbfs, nvtbf)
        push!(bfTs, vbfT)

        println("Done!")
        println("共得到 $ntetra 个四面体， $nvtbf 个 $vbfT 基函数。")

    end #begin

    # 有六则创建六、基函数信息
    meshData.hexanum   > 0     && begin
        # 开始构造六信息、基函数信息
        print("Constructing VIE basis function based on hexahedras...\t")
        hexameshData   =   HexahedraMesh( meshData.hexanum,    meshData.node,  meshData.hexahedras)
        if (vbfT === :PWC) || (vbfT === :RBF)
            hexasInfo, vhbfsInfo =   getHexasInfo(hexameshData, vbfT)
        else
            throw("体积分基函数类型出错或与网格类型不匹配，请检查！")
        end
        # 设置介质体介电参数
        setGeosPermittivity!(hexasInfo)
        # 三角形数量
        nhexa   =   length(hexasInfo)
        # 基函数数量
        nvhbf   =   length(vhbfsInfo)
        # 混合网格时更新 id 信息
        @threads for geo in hexasInfo
            geo.hexaID  +=  ngeo
            geo.inBfsID .+=  nbf
        end
        # 混合网格时更新基函数信息
        @threads for bf in vhbfsInfo
            bf.bfID      +=  nbf
            for ii in 1:length(bf.inGeo)
                bf.inGeo[ii] == 0 && continue
                if vbfT === :RBF
                    bf.inGeo[ii] +=  ngeo
                elseif vbfT === :PWC
                    bf.inGeo    +=  ngeo
                end
            end
        end
        # 几何信息更新
        ngeo    +=   nhexa
        nbf     +=   nvhbf

        push!(nbfs, nvhbf)
        push!(bfTs, vbfT)

        println("Done!")
        println("共得到 $nhexa 个六面体， $nvhbf 个 $vbfT 基函数。")
        
    end #begin

    # 更新保存的参数信息
    open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
        if all(isequal(:PWC), bfTs)
            record_BFsInfo(:PWC, nbf; io = f)
        else
            record_BFsInfo(bfTs, nbfs; io = f)
        end
    end

    ## 根据情况返回相关值
    if (meshData.trinum  ==  meshData.geonum)
        return ntri, nsbf, trisInfo, stbfsInfo
    elseif (meshData.tetranum  ==  meshData.geonum)
        return ntetra, nvtbf, tetrasInfo, vtbfsInfo
    elseif (meshData.hexanum  ==  meshData.geonum)
        return nhexa, nvhbf, hexasInfo, vhbfsInfo
    else
        geosInfo = Vector[]
        bfsInfo  = Vector[]
        if (meshData.trinum > 0)
            push!(geosInfo, [trisInfo])
            push!(bfsInfo, [stbfsInfo])
        end
        if (meshData.tetranum > 0)
            push!(geosInfo, [OffsetVector(tetrasInfo, meshData.trinum)])
            push!(bfsInfo, [vtbfsInfo])
        end
        if (meshData.hexanum > 0)
            push!(geosInfo, [OffsetVector(hexasInfo, meshData.trinum + meshData.tetranum)])
            push!(bfsInfo, [vhbfsInfo])
        end
        return ngeo, nbf, geosInfo, bfsInfo
    end
end

"""
通过文件名直接读取网格元
输入:
meshFileName::  文件名字符串
meshUnit::      网格单位
sbfT    ::      面基函数类型
vbfT    ::      体基函数类型
"""
function getCellsFromFileName(meshFileName; meshUnit = MeshUnit, sbfT = :RWG, vbfT = :nothing)
    # 更新仿真参数
    modiSimulationParams!(;meshfilename = meshFileName, meshunit = meshUnit)

    updateVSBFTypes!(;sbfType = sbfType, vbfType = vbfType)
    # 保存仿真参数
    saveSimulationParams(;meshfilename = meshFileName, sbfT = sbfT, vbfT = vbfT)

    # 运行函数得到网格数据
    meshData    =   getMeshData(meshFileName,  meshUnit = meshUnit)

    println("\tDone!")

    return meshData

end

"""
通过文件名直接读取网格元、创建基函数信息
输入:
meshFileName::  文件名字符串
meshUnit::      网格单位
sbfT    ::      面基函数类型
vbfT    ::      体基函数类型
"""
function getCellsBFsFromFileName(meshFileName; meshUnit = MeshUnit, sbfT = :RWG, vbfT = :nothing)
    
    # 更新仿真参数
    modiSimulationParams!(;meshfilename = meshFileName, meshunit = meshUnit)
    # 运行函数得到网格数据
    meshData    =   getCellsFromFileName(meshFileName; meshUnit = meshUnit, sbfT = sbfT, vbfT = vbfT)

    println("\tDone!")

    # 计算网格元、基函数
    ngeo, nbf, geosInfo, bfsInfo = getCellsBFs(meshData, vbfT)

    ngeo, nbf, geosInfo, bfsInfo
end

"""
通过文件名直接读取网格元
输入:
meshFileName::  文件名字符串
meshUnit::      网格单位
sbfT    ::      面基函数类型
vbfT    ::      体基函数类型
"""
function getBFsFromMeshData(meshData; sbfT = :nothing, vbfT = :nothing)
    # 保存仿真参数
    saveSimulationParams(;sbfT = sbfT, vbfT = vbfT)
    # 更新面、体基函数类型
    updateVSBFTParams!(;sbfT = sbfT, vbfT = vbfT)

    # 计算网格元、基函数
    @clock "构建网格元、基函数" begin
        ngeo, nbf, geosInfo, bfsInfo = getCellsBFs(meshData, vbfT)
    end
    memory["网格元"]    =   Base.summarysize(geosInfo)
    memory["基函数"]    =   Base.summarysize(bfsInfo)
    return  ngeo, nbf, geosInfo, bfsInfo

end