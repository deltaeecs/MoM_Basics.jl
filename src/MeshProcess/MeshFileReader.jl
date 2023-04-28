module MeshFileReader

using FiniteMesh, Meshes


export  readMesh, 
        generateMeshFromFile, 
        fromNodesCells2SimpleMesh, 
        fromNodesConns2SimpleMesh


"""
读取网格文件并根据输入的单位将坐标单位换算成 "米" ，注意对于 Hypermesh 生成的网格，其单位多为毫米(mm)
输入：
fn:: 文件路径 + 名
unit:: 网格文件所采用的单位，目前支持 米(m)， 厘米(cm)， 毫米(mm) 三种单位" 
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

"""
给出三种 "FiniteMesh" 读取出的网格类型到对应的 "Meshes" 网格类型的字典
"""
const cellTypesDict = Dict("triangle" => Triangle, "tetra" => Tetrahedron, "hexahedron" => Hexahedron)
const cellTypes     = Dict(3 => Triangle, 4=> Tetrahedron, 8=> Hexahedron)

# 测试
# cells, nodes = read_mesh("meshfiles/cone_1GHz.nas", :mm)


"""
用读出的数据生成网格
"""
function generateMeshFromFile(fn::T, unit::Symbol) where{T<:AbstractString}
    # 读取单元和点
    cells, nodes = readMesh(fn, unit)
    # 由单元和点成网格
    fromNodesCells2SimpleMesh(nodes, cells)

end

"""
nodes and cells to mesh
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
    mesh = SimpleMesh(points, cons)

end # function


"""
nodes and connects to mesh
nodes:: Matrix{FT} (3, n)
connects:: Matrix{IT} (ngeo, ncell)
geoTi 1:triangle, 2:tetrahedra, 3:hexahedra
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