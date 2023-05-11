module MeshFileReader

using FiniteMesh, Meshes

export  readMesh, 
        generateMeshFromFile, 
        fromNodesCells2SimpleMesh, 
        fromNodesConns2SimpleMesh

"""
    readMesh(fn::T, unit::Symbol) where{T<:AbstractString}

读取名称为 `fn` 的网格文件并根据输入的单位将坐标单位换算成 `米(m)` ，
注意对于 Hypermesh 生成的网格，一般不存在单位，可根据实际情况填写单位。 
"""
function readMesh(fn::T, unit::Symbol) where{T<:AbstractString}
    
    cells, nodes::Matrix{Float64} = read_mesh(fn)

    if unit == :mm
        nodes .*= 1e-3
    elseif unit == :cm
        nodes .*= 1e-2
    elseif unit == :m
        nothing
    else
        throw("目前只支持 米(m)， 厘米(cm)， 毫米(mm) 三种单位")
    end

    cells, nodes

end

# """
# [`FiniteMesh`](@ref) 读取出的三角形、四面体、六面体网格类型到对应的 [`Meshes`](@ref) 网格类型的字典。
# """
# const cellTypesDict = Dict("triangle" => Triangle, "tetra" => Tetrahedron, "hexahedron" => Hexahedron)
# const cellTypes     = Dict(3 => Triangle, 4=> Tetrahedron, 8=> Hexahedron)

"""
    generateMeshFromFile(fn::T, unit::Symbol) where{T<:AbstractString}

用读出的数据生成 [`Meshes`](@ref) 网格。
"""
function generateMeshFromFile(fn::T, unit::Symbol) where{T<:AbstractString}
    # 读取单元和点
    cells, nodes = readMesh(fn, unit)
    # 由单元和点成网格
    fromNodesCells2SimpleMesh(nodes, cells)

end

"""
    fromNodesCells2SimpleMesh(nodes, cells)

用读出的网格 `cells` 和节点 `nodes` 数据生成 [`Meshes`](@ref) 网格。
"""
function fromNodesCells2SimpleMesh(nodes, cells)
    # 点数组到点类型
    points  =   Point3f.(eachrow(nodes))
    
    ## 复合类型网格有问题，
    "=
    # 连接
    # cons    =   Connectivity[]
    # 处理不同类型的网格
    # for ii in 1:length(cells.type)
    #     con =  connect.(Tuple.(eachrow(cells.index[ii])), cellTypes[cells.type[ii]])
    #     push!(cons, con...)
    # end
    ="
    # 将就处理单类型网格
    cons =  connect.(Tuple.(eachrow(cells.index[1])), cellTypesDict[cells.type[1]])

    # 构建网格
    return SimpleMesh(points, cons)

end # function


"""
    fromNodesConns2SimpleMesh(nodes, cons)

用节点 `nodes` 和连接 `cons` 数据生成 [`Meshes`](@ref) 网格。
"""
function fromNodesConns2SimpleMesh(nodes, cons)
    # 点数组到点类型
    points  =   Point3.(eachcol(nodes))

    ngeo    =   size(cons, 1)

    # 将就处理单类型网格
    consMesh =  connect.(Tuple.(eachcol(cons)), cellTypes[ngeo])

    # 构建网格
    mesh = SimpleMesh(points, consMesh)

end # function


end # module