# 本文件构建六面体形类和一些积分函数、奇异值处理等

"""
    HexahedraMesh{IT, FT} <: MeshDataType

六面体网格文件，包括六面体数 `hexnum`、节点 `node`、构成六面体的节点 id 数组 `hexahedras` 等信息。
"""
struct HexahedraMesh{IT, FT} <: MeshDataType
    hexnum      ::Int
    node        ::Array{FT, 2}
    hexahedras  ::Array{IT, 2}
end

include("Quadrangle.jl")

# """
# 重载 getproperty 函数用于直接提取六面体信息
# """
# function Base.getproperty(tetra::HexahedraInfo, sym::Symbol)
#     if sym === :faces
#         return view(uniqueTris, tetra.facesID)
#     else # fallback to getfield
#         return getfield(tetra, sym)
#     end
# end

"""
    HexahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex} <: VolumeCellType{IT, FT, CT}

六面体信息结构体：
```
hexaID      ::IT                    编号
volume      ::FT                    体积
ε           ::CT                    相对介电常数
κ           ::CT                    介质对比度
center      ::MVec3D{FT}            中心坐标
verticesID  ::MVector{8, IT}        所在节点id
vertices    ::MMatrix{3, 8, FT,24}  六面体4个角点坐标，每列为一个点
facesn̂      ::MMatrix{3, 8, FT,24}  四个面的外法向量
facesvid    ::MMatrix{3, 8, IT,24}  四个面包含的四个id
facesArea   ::MVector{6, FT}        四个面的面积（根据为unitri正负部分赋予正负号）
faces       ::Vector{Quads4Hexa{IT, FT}}    四个面的具体信息
inBfsID     ::Vector{IT}            六面体所在的基函数的ID
```
"""
mutable struct HexahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex} <: VolumeCellType{IT, FT, CT}
    hexaID      ::IT
    volume      ::FT
    ε           ::CT
    κ           ::CT
    center      ::MVec3D{FT}
    verticesID  ::MVector{8, IT}
    vertices    ::MMatrix{3, 8, FT, 24}
    facesn̂      ::MMatrix{3, 6, FT, 18}
    facesArea   ::MVector{6, FT}
    faces       ::Vector{Quads4Hexa{IT, FT}}
    inBfsID     ::Vector{IT}
end

"""
    HexahedraInfo{IT, FT, CT}(hexaID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat, CT<:Complex}

`HexahedraInfo` 的默认构造函数，除了输入的编号 `hexaID` 外所有元素置零。
"""
function HexahedraInfo{IT, FT, CT}(hexaID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat, CT<:Complex}
    volume::FT  =   zero(FT)
    ε           =   one(CT)
    κ           =   zero(CT)
    center      =   zero(MVec3D{FT})
    verticesID  =   zero(MVector{8, IT})
    vertices    =   zero(MMatrix{3, 8, FT, 24})
    facesn̂      =   zero(MMatrix{3, 6, FT, 18})
    facesArea   =   zero(MVector{6, FT})
    faces       =   Vector{Quads4Hexa{IT, FT}}(undef, 6)
    inBfsID     =   zeros(IT, 6)
    return HexahedraInfo{IT, FT, CT}(hexaID, volume, ε, κ, center, verticesID, vertices,
                                        facesn̂, facesArea, faces, inBfsID)
end


## 修改的话与四边形一起修改，保证与四边形的对齐（同一纬度采样点一样）
```
正常求积六面体高斯求积点数。
```
const GQPNHexa      =   8
```
处理奇异性时的六面体高斯求积点数。
```
const GQPNHexaSglr  =   27
```
处理超高奇异性时的六面体高斯求积点数。
```
const GQPNHexaSSglr =   64

```
正常求积六面体高斯求积信息。
```
const HexaGQInfo      =   GaussQuadrature4Geo.GaussQuadratureInfo(:Hexahedron, GQPNHexa, Precision.FT)
```
处理奇异性时六面体高斯求积信息。
```
const HexaGQInfoSglr  =   GaussQuadrature4Geo.GaussQuadratureInfo(:Hexahedron, GQPNHexaSglr, Precision.FT)
```
处理超奇异性时六面体高斯求积信息。
```
const HexaGQInfoSSglr =   GaussQuadrature4Geo.GaussQuadratureInfo(:Hexahedron, GQPNHexaSSglr, Precision.FT)


"""
    setHexaCoor!( hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}, hexaMeshData::HexahedraMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

在预分配好的六面体数组 `hexasInfo` 里写入 `hexaMeshData` 中对应的六面体编号、点坐标、中心位置数据。
"""
function setHexaCoor!( hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}, 
                        hexaMeshData::HexahedraMesh{IT, FT}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}
    @threads for i in 1:length(hexasInfo)
        hexasInfo[i].hexaID      =  i
        hexasInfo[i].verticesID  .=  view(hexaMeshData.hexahedras, :, i)
        for ii in 1:8
            hexasInfo[i].vertices[:,ii]   .=  @view hexaMeshData.node[:,hexasInfo[i].verticesID[ii]]
        end
        hexasInfo[i].center   =   mean(hexasInfo[i].vertices, dims = 2)
    end

    return nothing
end


"""
    setHexaParam!(hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

计算六面体体积、面外法向量、面积，并写入 `hexasInfo` 。
"""
function setHexaParam!(hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}

    @threads for i in 1:length(hexasInfo)
        @inbounds begin 
        # 第 i 六面体
        hexaInfo    =   hexasInfo[i]
        # 六面体的网格
        hexa    =   Hexahedron(Meshes.Point3.(eachcol(hexaInfo.vertices)))
        # 体积
        hexaInfo.volume =   measure(hexa)
        # 六个四边形面
        quads   =   boundaryRBF(hexa)
        # 对四个面循环计算面的外法向量，注意第 facei 个面为 第 facei 个 vertice 对面的四个点构成的四边形
        for facei in 1:6
            #  第 facei 个四边形
            quadi   =  quads[facei]
            # 四个点的视图
            faceVertViews   =   quadi.vertices

            # 计算面法向量（未定向\未归一化）
            @views facen̂    =   hexaInfo.facesn̂[:, facei]
            @views facen̂   .=   cross(faceVertViews[2] - faceVertViews[1],  faceVertViews[2] - faceVertViews[3])
            # 未归一化面法向量的模即为四边形面积两倍
            faceArea    =   area(quadi)
            # 归一化面法向量
            facen̂     ./=   norm(facen̂)

            # 定向，第 facei 个 vertice 点到其余四个点中的一个的向量与外法向量的点乘结果应该大于0
            @views ((faceVertViews[1].coords .- hexaInfo.center) ⋅ facen̂ < 0) && begin facen̂ .*= -1 end

            # 写入结果
            hexaInfo.facesArea[facei]  =   faceArea
            # 面法向量结果采用的是内存视图，不必再重新写入
        end # for facei
    end # begin
    end # for i
    return nothing

end # function


"""
    getHexasInfo(hexameshData::HexahedraMesh{IT, FT}, VolumeBFType::Symbol) where{IT, FT}

根据六面体网格信息 `hexameshData` 和体基函数类型 `VolumeBFType` 生成网格信息向量 `hexasInfo` 和基函数信息向量 `bfsInfo` 。
"""
function getHexasInfo(hexameshData::HexahedraMesh{IT, FT}, VolumeBFType::Symbol) where{IT, FT}

    # # 先生成 SimpleMesh
    # hexaSimpleMesh = MeshFileReader.fromNodesConns2SimpleMesh(hexameshData.node, hexameshData.hexahedras)
    # 分配四边形信息类型数组
    hexasInfo       =   [HexahedraInfo{IT, FT, Complex{FT}}()  for _ in 1:hexameshData.hexnum ]
    # 写入六面体的八个点的id和点坐标、几何中心数据
    setHexaCoor!(hexasInfo, hexameshData)
    # 写入几何中心、边界的四边形面积、外法向量、面法向量、面积
    setHexaParam!(hexasInfo)
    # 计算基函数信息并修改四边形信息中基函数相关项
    bfsInfo  =   setQuad4Hexas!(hexameshData, hexasInfo, Val(VolumeBFType))
    return hexasInfo, bfsInfo
end


"""
    setδκ!(hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}

设置六面体网格信息 `hexasInfo` 中每个面上的介质对比度差值。
"""
function setδκ!(hexasInfo::AbstractVector{HexahedraInfo{IT, FT, CT}}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    # 循环置零
    @threads for it in eachindex(hexasInfo)
        # 六面体
        hexa   =   hexasInfo[it]
        # 对面循环
        for iface in 1:6
            # 面
            face =  hexa.faces[iface]
            # 置零
            face.δκ = 0
        end
    end
    # 循环设置δκ
    @threads for it in eachindex(hexasInfo)
        # 六面体
        hexa   =   hexasInfo[it]
        # 六面体的介质对比度
        κ   =    hexa.κ
        # 对面循环
        for iface in 1:6
            # 面
            face =  hexa.faces[iface]
            # 根据在正负面决定加上或减去 κ
            temp =  (hexa.facesArea[iface] > 0 ? κ : -κ)
            # 写入值
            face.δκ += temp
        end
    end
end

@doc """
    getGQPHexa(hexa::HexahedraInfo, ii::IT) where {IT <: Integer}
    getGQPHexa(hexa::HexahedraInfo)

计算 `hexa` 正常求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPHexa(hexa::HexahedraInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds hexa.vertices * HexaGQInfo.coordinate[:, i]
end
function getGQPHexa(hexa::HexahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    hexa.vertices * HexaGQInfo.coordinate
end

@doc """
    getGQPHexaSglr(hexa::HexahedraInfo, ii::IT) where {IT <: Integer}
    getGQPHexaSglr(hexa::HexahedraInfo)

计算 `hexa` 处理奇异性求积的第 `i` 个或所有高斯求积坐标。
"""
function getGQPHexaSglr(hexa::HexahedraInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds hexa.vertices * HexaGQInfoSglr.coordinate[:, i]
end
@inline function getGQPHexaSglr(hexa::HexahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    hexa.vertices * HexaGQInfoSglr.coordinate
end

@doc """
    getGQPHexaSSglr(hexa::HexahedraInfo, ii::IT) where {IT <: Integer}
    getGQPHexaSSglr(hexa::HexahedraInfo)

计算 `hexa` 处理超奇异性求积的第 `i` 个或所有高斯求积坐标。
"""
@inline function getGQPHexaSSglr(hexa::HexahedraInfo, i::IT) where {IT <: Integer}
    # 直接向量相乘，此处采用的是矩阵乘法
    @views @inbounds hexa.vertices * HexaGQInfoSSglr.coordinate[:, i]
end
@inline function getGQPHexaSSglr(hexa::HexahedraInfo)
    # 直接向量相乘，此处采用的是矩阵乘法
    hexa.vertices * HexaGQInfoSSglr.coordinate
end