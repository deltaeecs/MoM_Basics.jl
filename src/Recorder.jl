"""
程序计时器字典
"""
const timer = Dict{String, Float64}()

"""
程序内存记录字典
"""
const memory= Dict{String, Float64}()


"""
    initialize_time_and_mem()

初始化计时器 `timer` 和内存记录 `memory`。
"""
function initialize_time_and_mem()

    map(k -> delete!(timer, k), collect(keys(timer)))
    map(k -> delete!(memory, k), collect(keys(memory)))

    return nothing

end

"""
    clock(message, ex)
	
将表达式 `ex` 的运行时间以 `message` 为键保存在字典 `timer` 中。
"""
macro clock(message, ex)
    quote
        local t_0 = time()
        $(esc(ex))
        local t_1 = time()
        
        timer[$message] = t_1 - t_0
    end
end

"""
    byte2other(v, mem_unit)

将以字节为单位的内存数据 `v` 转换为其它单位 `mem_unit`。
"""
function byte2other(v, mem_unit)
    return if mem_unit == :B
        v
    elseif mem_unit == :KB
        v/1024
    elseif mem_unit == :MB
        v/1024/1024
    elseif mem_unit == :GB
        v/1024/1024/1024
    elseif mem_unit == :TB
        v/1024/1024/1024/1024
    end
end

"""
    second_to_other(v, time_unit)

将以秒为单位的时间数据 `v` 转换为其它单位 `time_unit`。
"""
function second_to_other(v, time_unit)
    return if time_unit == :s
        v
    elseif time_unit == :ms
        v*1e3
    elseif time_unit == :us
        v*1e6
    elseif time_unit == :ns
        v*1e9
    end
end

"""
    show_memory_time(io::IO=Core.stdout; mem_unit = :MB, time_unit = :s)

展示时间内存消耗数据记录。
"""
function show_memory_time(io::IO=Core.stdout; mem_unit = :MB, time_unit = :s)
    @printf io "%19s\n" "内存"
    for (k, v) in memory
        m = byte2other(v, mem_unit)
        @printf(io,  "%-20s %10.3f %2s\n", k, m, mem_unit)
    end
    @printf io "%19s\n" "时间"
    for (k, v) in timer
        t = second_to_other(v, time_unit)
        @printf(io,  "%-20s %10.3f %2s\n", k, t, time_unit)
    end
    nothing
end


@doc """
    record_CellInfo(io::IO = Core.stdout; ntri = 0, ntetra = 0, nhexa = 0)
    record_CellInfo(meshData; io::IO = Core.stdout)

在 `io` 中记录网格单元数量 `ntri, ntetra, nhexa`。
"""
function record_CellInfo(io::IO = Core.stdout; ntri = 0, ntetra = 0, nhexa = 0)

    @printf io "%21s\n" "网格信息"
    ntri > 0 && @printf io "%-20s %13i\n" "Triangles" ntri
    ntetra > 0 && @printf io "%-20s %13i\n" "Tetrahedrons" ntetra
    nhexa > 0 && @printf io "%-20s %13i\n" "Hexahedrons" nhexa

    nothing

end
function record_CellInfo(meshData; io::IO = Core.stdout)
    record_CellInfo(io; ntri = meshData.trinum, ntetra = meshData.tetranum, nhexa = meshData.hexanum)
end

@doc """
    record_BFsInfo(bfT::Symbol, nbf::Int; io::IO = Core.stdout)

在 `io` 中记录基函数类型 `bfT` 和数量 `nbf`。
"""
function record_BFsInfo(bfT::Symbol, nbf::Int; io::IO = Core.stdout)

    @printf io "%22s\n" "基函数信息"
    nbf > 0 && @printf io "%-20s %13i\n" bfT nbf

    nothing

end
function record_BFsInfo(bfTs::Vector{Symbol}, nbfs::Vector{Int}; io::IO = Core.stdout)

    @printf io "%22s\n" "基函数信息"
    @printf io "%-21s %13i\n" "总的基函数数量" sum(nbfs)
    for (bfT, nbf) in zip(bfTs, nbfs)
        bfT != :nothing && @printf io "%-20s %13i\n" bfT nbf
    end

    nothing

end

