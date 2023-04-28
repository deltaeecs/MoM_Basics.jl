#=本文件用于处理网格文件, 目前仅支持.nas格式的网格=#

"""构建网格文件抽象类型"""
abstract type MeshFormat end
"nas文件符合类型，存在字段pathname保存文件名称和相对路径"
struct NasMesh<:MeshFormat
    pathname ::String

    function NasMesh(MeshFile::String)
        # nas格式
        !occursin(".nas",  MeshFile) &&  throw("输入的不是 .nas 文件，目前只支持该类型网格文件。")
        new(MeshFile)
    end
end

abstract type MeshDataType end

"""
网格文件类型，目前包含三角形\\四面体\\六面体信息，包含属性如下：
geonum      :: 包含的所有网格元的数量
meshT       :: 网格类型 单一的以其网格类型表示，混合以 VSCellType 表示
trinum      : 包含的三角形数量
tetranum    : 包含的四面体数量
hexanum     : 包含的六面体数量
node        : 节点坐标数组:(3*nodenum)
triangles   : 三角形包含的nodeid数组:(3*trinum)
tetrahedras : 四面体包含的nodeid数组:(4*tetranum)
hexahedras  : 六面体包含的nodeid数组:(6*hexanum)
"""
struct MeshNodeTriTetraHexa{IT, FT} <: MeshDataType
    geonum      ::Int
    meshT       ::DataType
    trinum      ::Int
    tetranum    ::Int
    hexanum     ::Int
    node        ::Array{FT, 2}
    triangles   ::Array{IT, 2}
    tetrahedras ::Array{IT, 2}
    hexahedras  ::Array{IT, 2}
end

"""
读取文件中的节点坐标、三角形点、四面体点
"""
function getNodeTriTetraFekoNas(pathname::ST, FT::Type{T}=Precision.FT) where {ST <: AbstractString,T<:AbstractFloat}
    # 更新仿真参数
    modiSimulationParams!(;meshfilename = pathname)
    meshfile    = NasMesh(pathname)
    # 节点数、三角形数、四面体数
    local nodenum     = zero(Int64)
    local trinum      = zero(Int64)
    local tetranum    = zero(Int64)
    # 打开文件, 找出点数
    open(meshfile.pathname, "r") do NasMesh
        for line in eachline(NasMesh)
            contents  =   split(line)
            occursin("GRID*",       line) && (nodenum += 1;                         continue)
            occursin("triangles",   line) && (trinum   = parse(Int32, contents[3]); continue)
            occursin("tetrahedra",  line) && (tetranum = parse(Int32, contents[3]); continue)
        end
    end
    # 精度
    FloatDtype  =   Precision.FT
    # 预分配存储节点、三角形、四面体数组
    # 无符号的索引速度比有符号慢不少，故不次采用
    # LocalInt64=   ((trinum > typemax(Int32) | tetranum > typemax(Int32)) ?  Int64 : Int32)
    node        =   zeros(FloatDtype,   (3, nodenum ))
    triangles   =   zeros(Int64,     (3, trinum  ))
    tetrahedras =   zeros(Int64,     (4, tetranum))

    # 重新打开文件读入以上几个信息
    open(meshfile.pathname, "r") do NasMesh
        NodeID = 1; TriID = 1; TetraID = 1;
        for line in eachline(NasMesh)
            contents = split(line)
            if      contents[1] == "GRID*"
                node[1:2, NodeID]    .=   try
                        [parse(FloatDtype, contents[3]), parse(FloatDtype, contents[4])]
                    catch ArgumentError
                        length(contents[3]) == 32 ? [parse(FloatDtype, contents[3][1:16]), parse(FloatDtype, contents[3][17:end])] : [parse(FloatDtype, contents[3][1:15]), parse(FloatDtype, contents[3][16:end])]
                    end
            elseif  contents[1] == "*"
                node[3,   NodeID]    =    try 
                    parse(FloatDtype, contents[3])
                catch ArgumentError
                    -parse(FloatDtype, split(contents[2], "-"; limit = 2)[2])
                end
                NodeID  +=  1
            elseif  contents[1] == "CTRIA3"
                triangles[:, parse(Int64, contents[2])]  .=  
                [parse(Int64, contents[4]), parse(Int64, contents[5]), parse(Int64, contents[6])]
                TriID   +=  1
            elseif  contents[1] == "CTETRA"
                tetrahedras[:, parse(Int64, contents[2])].=   [
                parse(Int64, contents[4]), parse(Int64, contents[5]), 
                parse(Int64, contents[6]), parse(Int64, contents[7])]
            end
        end
    end
    # 总的网格单元数
    geonum  =   trinum + tetranum
    # 网格文件类型
    meshT   =   if (trinum > 0) && (tetranum > 0)
        VSCellType
    elseif (trinum > 0)
        TriangleInfo
    elseif (tetranum > 0)
        TetrahedraInfo
    else
        throw("输入网格文件出错，请检查！")
    end
    # 返回复合类型MeshNodeTriTetra
    MeshNodeTriTetraHex{Int64, FloatDtype}(geonum, meshT, trinum, tetranum, 0, node, triangles, tetrahedras, zeros(0, 0))
end

"""
读取文件中的节点坐标、三角形点、四面体点
"""
function getNodeTriTetraHexaNas(pathname::ST; FT::Type{T}=Precision.FT, meshUnit = :mm) where {ST <: AbstractString,T<:AbstractFloat}
    # 更新仿真参数
    modiSimulationParams!(;meshfilename = pathname, meshunit = meshUnit)
    # Nas网格
    meshfile    = NasMesh(pathname)

    # 打开文件, 找出# 节点数、三角形数、四面体数、六面体数
    linenum, nodenum, trinum, tetranum, hexanum = open(meshfile.pathname, "r") do NasMesh
        linenum = 0; nodenum = 0; trinum  = 0;tetranum= 0;hexanum = 0
        while true
            linenum += 1
            line = readline(NasMesh)
            startswith(line, r"\$|//|#") && continue
            if occursin("GRID",    line)
                (nodenum   += 1; continue)
            elseif occursin("CTRIA3",  line)
                (trinum    += 1; continue)
            elseif occursin("CTETRA",  line)
                (tetranum  += 1; continue)
            elseif occursin("CHEXA",   line)
                (hexanum   += 1; continue)
            else
                eof(NasMesh) ? break : continue
            end
        end
        linenum, nodenum, trinum, tetranum, hexanum
    end
    

    # 预分配存储节点、三角形、四面体、六面体数组
    # 无符号的索引速度比有符号慢不少，故不次采用
    # LocalInt64=   ((trinum > typemax(Int32) | tetranum > typemax(Int32)) ?  Int64 : Int32)
    # 点的 id 不一定按顺序来，因此早建立一个新索引
    nodeO2LID   =   Dict{Int, Int}()
    node        =   zeros(FT,(3, nodenum ))
    triangles   =   zeros(Int64,     (3, trinum  ))
    tetrahedras =   zeros(Int64,     (4, tetranum))
    hexahedras  =   zeros(Int64,     (8, hexanum ))

    # 重新打开文件读入以上几个信息
    open(meshfile.pathname, "r") do f
        NodeID = 0; TriID = 0; TetraID = 0; HexaID = 0
        nline = readline(f)
        # 进度条
        pmeter = Progress(linenum, "处理网格文件中...")
        contents = String[]
        while true
            next!(pmeter)
            line        =   nline
            nline       =   readline(f)
            # 空行或注释
            if (!eof(f)  && length(line) < 9) || startswith(line, r"\$|//|#")
                contents = String[]
                continue
            end

            # 连续行或新行
            if startswith(line, r"\+|\*")
                append!(contents, _chunk_line(line)[2:end])
            else
                contents = _chunk_line(line)
            end
            # 下一行连续
            startswith(nline, r"\+|\*") && continue

            # 开始解析
            emptyc = isempty(contents)
            if  !emptyc && startswith(contents[1], "GRID")
                NodeID += 1
                nodeO2LID[parse(Int, contents[2])] = NodeID
                node[1, NodeID] =   _nastran_string_to_float(contents[4])
                node[2, NodeID] =   _nastran_string_to_float(contents[5])
                node[3, NodeID] =   _nastran_string_to_float(contents[6])

            elseif  !emptyc && startswith(contents[1], "CTRIA3")
                TriID += 1
                triangles[1, TriID]     =   nodeO2LID[parse(Int64, contents[4])]
                triangles[2, TriID]     =   nodeO2LID[parse(Int64, contents[5])]
                triangles[3, TriID]     =   nodeO2LID[parse(Int64, contents[6])]
            elseif  !emptyc && startswith(contents[1], "CTETRA")
                TetraID += 1
                tetrahedras[1, TetraID] =   nodeO2LID[parse(Int64, contents[4])]
                tetrahedras[2, TetraID] =   nodeO2LID[parse(Int64, contents[5])]
                tetrahedras[3, TetraID] =   nodeO2LID[parse(Int64, contents[6])]
                tetrahedras[4, TetraID] =   nodeO2LID[parse(Int64, contents[7])]
            elseif !emptyc && startswith(contents[1], "CHEXA")
                HexaID += 1
                hexahedras[1, HexaID]   =   nodeO2LID[parse(Int64, contents[4])]
                hexahedras[2, HexaID]   =   nodeO2LID[parse(Int64, contents[5])]
                hexahedras[3, HexaID]   =   nodeO2LID[parse(Int64, contents[6])]
                hexahedras[4, HexaID]   =   nodeO2LID[parse(Int64, contents[7])]
                hexahedras[5, HexaID]   =   nodeO2LID[parse(Int64, contents[8])]
                hexahedras[6, HexaID]   =   nodeO2LID[parse(Int64, contents[9])]
                hexahedras[7, HexaID]   =   nodeO2LID[parse(Int64, contents[10])]
                hexahedras[8, HexaID]   =   nodeO2LID[parse(Int64, contents[11])]

            else
                # 其它情况不处理
                eof(f) ? break : continue
            end

            contents = String[]
        end # while
        # 
        println("网格文件处理完毕，共得到 $nodenum 个节点、$trinum 个三角形、$tetranum 个四面体、$hexanum 个六面体")
    end # open

    if meshUnit == :mm
        node .*= 1e-3
    elseif meshUnit == :cm
        node .*= 1e-2
    elseif meshUnit == :m
        nothing
    else
        throw("目前只支持 米(m)， 厘米(cm)， 毫米(mm) 三种单位")
    end

    # 总的网格单元数
    geonum  =   trinum + tetranum + hexanum
    # 网格文件类型
    meshT   =   if (trinum > 0) && ((tetranum > 0) || (hexanum > 0))
        VSCellType
    elseif (trinum > 0)
        TriangleInfo{Int, FT}
    elseif (tetranum > 0)
        TetrahedraInfo{Int, FT, Complex{FT}}
    elseif (hexanum > 0)
        HexahedraInfo{Int, FT, Complex{FT}}
    else
        throw("输入网格文件出错，请检查！")
    end
    # 返回复合类型MeshNodeTriTetra
    MeshNodeTriTetraHexa{Int64, FT}(geonum, meshT, trinum, tetranum, hexanum, node, triangles, tetrahedras, hexahedras)

end

function _nastran_string_to_float(string)
    try
        return parse(Float64, string)
    catch ValueError
        string = strip(string)
        return parse(Float64, string[1] * replace(string[2:end], "+"=>"e+", "-"=>"e-"))
    end
end

function _chunk_line(line)
    length(line) < 8 && return strip.([line])
    chunks = begin
        if occursin(",", line)  # free format
            split(line, ",")
        else  # fixed format
            signs = strip(line[1:8])
            CHUNK_SIZE = if endswith(signs, '*') || startswith(signs, '*'); 16; else; 8; end
            [line[1:8], [line[i : CHUNK_SIZE + i - 1] for i in 9:CHUNK_SIZE:(length(line) - CHUNK_SIZE + 1)]...]
        end
    end # begin
    # everything after the 9th chunk is ignored
    re = strip.(chunks)
    (re[end] == "+" || re[end] == "*") && pop!(re)
    return re
end


"""
读取文件中的节点坐标、三角形点、四面体点、六面体点
"""
function getMeshData(meshFileName::String; meshUnit=:mm)
    # 创建结果目录
    ~isdir(SimulationParams.resultDir) && mkpath(SimulationParams.resultDir)
    if endswith(meshFileName, ".nas")
        @clock "网格文件读取" begin
            meshData = getNodeTriTetraHexaNas(meshFileName; meshUnit = meshUnit)
        end
        memory["网格文件"] = Base.summarysize(meshData)
        # 更新保存的参数信息
        open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
            record_CellInfo(meshData; io = f)
        end
        return meshData, nothing
    elseif endswith(meshFileName, ".dat")
        @clock "网格文件读取" begin
            meshData, εᵣs = getdatNodeElementParam(meshFileName, meshUnit=meshUnit)
        end
        memory["网格文件"] = Base.summarysize(meshData)
        # 更新保存的参数信息
        open(SimulationParams.resultDir*"/InputArgs.txt", "a+")  do f
            record_CellInfo(meshData; io = f)
        end
        return meshData, εᵣs
    end
end

"""
节点与网格元之间的连接稀疏矩阵。
"""
function getConnectionMatrix(meshData)
    # 节点数
    nNode   =   size(meshData.node, 2)
    # 网格元数
    nCell   =   meshData.geonum
    #  构建网格包含的节点 id 在稀疏矩阵中的行索引
    rowIndices  =   begin
        vcat(   reshape(meshData.triangles, :),
                reshape(meshData.tetrahedras, :),
                reshape(meshData.hexahedras, :))
    end
    # 构建网格包含的节点 id 在稀疏矩阵中的列指针
    colPtr      =   zeros(Int, nCell + 1)
    colPtr[1:end-1] .=   begin
        vcat(   1:3:3meshData.trinum,
                (1:4:4meshData.tetranum) .+ 3meshData.trinum, 
                (1:8:8meshData.hexanum) .+ (3meshData.trinum + 4meshData.tetranum))
    end
    colPtr[end] = length(rowIndices) + 1
    # 结果矩阵
    SparseMatrixCSC{Bool, Int}(nNode, nCell, colPtr, rowIndices, ones(Bool, length(rowIndices)))

end


"""
读取项目文件（.dat格式）
"""
function getdatNodeElementParam(pathname::ST; FT::Type{T}=Precision.FT, meshUnit = :m) where {ST <: AbstractString,T<:AbstractFloat}
    # 更新仿真参数
    modiSimulationParams!(;meshfilename = pathname, meshunit = meshUnit)
    
    # 打开文件
    mesh = open(pathname, "r")
    # 找出# 节点数、四面体数
    trinum  =   0
    hexanum =   0
    # 空过头两行
    for _ in 1:2 readline(mesh); end
    contents    =   split(readline(mesh))
    nodenum     =   parse(Int, contents[3])
    tetranum    =   parse(Int, contents[5])
    # 空过下一行
    for _ in 1:1 readline(mesh); end

    # 预分配存储节点、三角形、四面体、六面体数组
    # 无符号的索引速度比有符号慢不少，故不次采用
    # LocalInt64=   ((trinum > typemax(Int32) | tetranum > typemax(Int32)) ?  Int64 : Int32)
    # 点的 id 不一定按顺序来，因此早建立一个新索引
    node        =   zeros(FT,(3, nodenum ))
    triangles   =   zeros(Int64,     (3, trinum  ))
    tetrahedras =   zeros(Int64,     (4, tetranum))
    hexahedras  =   zeros(Int64,     (8, hexanum ))
    εᵣsNodes    =   zeros(Complex{FT},   nodenum )
    εᵣsTetra    =   zeros(Complex{FT},   tetranum )

    # 读取以上信息
    for i in 1:nodenum
        contents    =   split(readline(mesh))
        # 赋值
        for ixyz in 1:3
            node[ixyz, i] =   parse(FT, contents[ixyz])
        end
        εᵣsNodes[i] =   parse(FT, contents[5]) + parse(FT, contents[6])*1im
    end

    for i in 1:tetranum
        contents    =   split(readline(mesh))
        # 赋值
        for i_tetra_node in 1:4
            tetrahedras[i_tetra_node, i] = parse(Int, contents[i_tetra_node])
        end
        εᵣsTetra[i] = mean(i_node -> εᵣsNodes[i_node], view(tetrahedras, :, i))
    end

    close(mesh)

    if meshUnit == :mm
        node .*= 1e-3
    elseif meshUnit == :cm
        node .*= 1e-2
    elseif meshUnit == :m
        nothing
    else
        throw("目前只支持 米(m)， 厘米(cm)， 毫米(mm) 三种单位")
    end

    # 总的网格单元数
    geonum  =   trinum + tetranum + hexanum
    # 网格文件类型
    meshT   =   if (trinum > 0) && ((tetranum > 0) || (hexanum > 0))
        VSCellType
    elseif (trinum > 0)
        TriangleInfo{Int, FT}
    elseif (tetranum > 0)
        TetrahedraInfo{Int, FT, Complex{FT}}
    elseif (hexanum > 0)
        HexahedraInfo{Int, FT, Complex{FT}}
    else
        throw("输入网格文件出错，请检查！")
    end
    # 返回复合类型MeshNodeTriTetra
    MeshNodeTriTetraHexa{Int64, FT}(geonum, meshT, trinum, tetranum, hexanum, node, triangles, tetrahedras, hexahedras), εᵣsTetra

end
