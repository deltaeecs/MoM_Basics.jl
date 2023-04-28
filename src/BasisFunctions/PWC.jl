
# 定义 PWC 基函数数据结构
"""
PWC 基函数复合类型参数解释：
isbd    :   是否为边界元即半基函数，布尔类型
bfID   :   基函数编号，整形
inGeo   :   基函数所在两个四面体编号（半基函数为1个，赋值成一样的两个），长度为2的向量数组
center  :   面所在的三个点的中心，目前用途为在 Octree 中分组
"""
mutable struct PWC{IT<:Integer , FT<:AbstractFloat} <: ConstBasisFunction
    bfID        ::IT
    inGeo       ::IT
    center      ::MVec3D{FT}
end

"""
PWCbfstruct的默认构造函数，所有元素置零
"""
function PWC{IT, FT}() where {IT <: Integer, FT<:AbstractFloat}
    bfID        =    zero(IT)
    inGeo       =    zero(IT)
    center      =    zero(MVec3D{FT})
    return PWC{IT, FT}(bfID, inGeo, center)
end

PWC()   =    PWC{IntDtype, Precision.FT}()


"""
计算构成四面体的所有三角形，并写入四面体包含的四个三角形的id，给分片常数(PWC)基函数在赋值
"""
function setTriangles4Tetras!(tetrameshData::TetrahedraMesh{IT, FT}, tetrasInfo::Vector{TetrahedraInfo{IT, FT, CT}}, ::Val{:PWC}) where {IT, FT, CT}
    # 四面体数
    tetranum  =   tetrameshData.tetranum
    
    # 创建数组，用于后续保存所有面的点和面的对点信息(4)，然后是面所在的四面体id(1)，面所在 unitri(独立三角形) 编号(1)、正负等信息(2)，共 8 列
    face3opp1tri2unitri2   =   zeros(IT, (tetranum*4, 8))
    # 利用内存视图提取其中的所有面的点和面的对点信息信息
    face3opp1           =   reshape(view(face3opp1tri2unitri2, :, 1:4), (:, 16))
    # 利用内存视图提取四面体信息的tri2和基函数信息的unitri2
    tetra2              =   @view face3opp1tri2unitri2[:, 5:6]
    unitri2             =   @view face3opp1tri2unitri2[:, 7:8]

    ## 创建四面体索引( id+localid )并根据排序更新
    # 此维度为四面体id
    tetra2[:, 1]  .=   repeat(1:tetranum,  outer = 4 )
    # 此维度为该面在所在四面体的局部id
    tetra2[:, 2]  .=   repeat(1:4,  inner = tetranum )
    
    # 将每个四面体的四个面的点及其对点存放在 face3opp1 数组中
    #* 注意不可改变下面的顺序因为跟RWG基函数对应的自由点在三角形中顺序有关 *#
    face3opp1[:, 1] .=   face3opp1[:, 8] .=   face3opp1[:, 11] .=   face3opp1[:, 14]  .=  tetrameshData.tetrahedras[2, :]
    face3opp1[:, 2] .=   face3opp1[:, 5] .=   face3opp1[:, 12] .=   face3opp1[:, 15]  .=  tetrameshData.tetrahedras[3, :]
    face3opp1[:, 3] .=   face3opp1[:, 6] .=   face3opp1[:,  9] .=   face3opp1[:, 16]  .=  tetrameshData.tetrahedras[4, :]
    face3opp1[:, 4] .=   face3opp1[:, 7] .=   face3opp1[:, 10] .=   face3opp1[:, 13]  .=  tetrameshData.tetrahedras[1, :]
    # 重塑形状将每个边点组合放入同一维度
    face3opp1        =   reshape(face3opp1, (:, 4))

    ##  对前三列排序使得同一面的表示形式相同
    #   对前三个构成边的点排序
    @views sort!(face3opp1tri2unitri2[:, 1:3], dims = 2)
    # 将每个四面体的四个角标点作为整体排序，次序按v1,v2,v3
    face3opp1tri2unitri2   .=  sortslices(face3opp1tri2unitri2, dims = 1, by = x -> (x[1], x[2], x[3]), alg = ThreadsX.MergeSort)

    ## for 循环计算RWG基函数id
    @inbounds let unitriplusid   =   1
        unitri2[1, 1]  =   1; unitri2[1, 2] = 1
        # 初始值即1的位置要预计算一下
        unitriplusid   +=  1

        # 从2开始循环，跟前一边不同为 unitri+，跟前一边不同为 unitri-
        for ii in 2:(tetranum*4)
            if  (face3opp1[ii, 1]  != face3opp1[ii-1, 1]) ||
                (face3opp1[ii, 2]  != face3opp1[ii-1, 2]) ||
                (face3opp1[ii, 3]  != face3opp1[ii-1, 3])
                unitri2[ii, 1]      =   unitriplusid
                unitri2[ii, 2]      =   1
                unitriplusid       +=   1
            elseif  (face3opp1[ii, 1]  ==   face3opp1[ii-1, 1]) &&
                (face3opp1[ii, 2]  ==   face3opp1[ii-1, 2]) && 
                (face3opp1[ii, 3]  ==   face3opp1[ii-1, 3])
                unitri2[ii, 1]      =   unitriplusid  - 1
                unitri2[ii, 2]      =   -1
            end #if
        end #for
    end #begin

    # unitri函数数量( 同时也是独立三角形的数量 )
    unitrinum   =   unitri2[end, 1]
    # 单独三角形所在的 face3opp1 的行 id
    unitri2nodeRow    =   zeros(IT, unitrinum)
    iunitri     =   1
    for inode in 1:(tetranum*4)
        (unitri2[inode, 2] == 1) && begin
            unitri2nodeRow[iunitri]     =   inode
            iunitri += 1
        end
    end
    # 预分配所有的构成四面体的三角形
    uniqueTris  =   [Tris4Tetra{IT, FT}()  for _ in 1:unitrinum ]

    # 设置 uniqueTris 的三个点、边长、边法向量等信息
    @threads for iTri in 1:unitrinum
        # 该三角形
        tri =   uniqueTris[iTri]
        vertices    = tri.vertices

        # 该三角形对应的三个角点id
        trivertIDs  =   view(face3opp1, unitri2nodeRow[iTri], :)
        # 设置节点
        for iv in 1:3
            vertices[:, iv]     .=   view(tetrameshData.node, :, trivertIDs[iv])
        end

        # 边指向向量
        edgev̂          =   tri.edgev̂
        @views for edgei in 1:3
            edgev̂[:, edgei]    =   vertices[:, EDGEVpINTriVsID[edgei]]  - vertices[:, EDGEVmINTriVsID[edgei]]
        end

        # 边长
        for edgei in 1:3
            @views tri.edgel[edgei]   =   norm(edgev̂[:, edgei])
        end

        @views facen̂    =  -cross(edgev̂[:, 1], edgev̂[:, 3])
        facen̂    ./=   norm(facen̂)
        # 计算边的指向单位向量、外法单位向量
        edgen̂                  =   tri.edgen̂
        for edgei in 1:3
            # 边的指向向量单位化
            @views edgev̂[:, edgei]  ./=   tri.edgel[edgei]
            # 最后计算边外法向量
            @views edgen̂[:, edgei]    =   cross(edgev̂[:, edgei], facen̂)
        end

    end

    # 写入公共面数信息
    @views let faceii = 1, faceID::IT
        while faceii < 4*tetranum
            faceID   =   unitri2[faceii, 1]
            # 挑出第2[faceii,1]个面信息
            uniqueTri   =   uniqueTris[faceID]
            if unitri2[faceii, 1] ==  unitri2[faceii + 1, 1]
                # 全基函数信息
                uniqueTri.isbd  =   false
                # 跳过负基函数部分
                faceii += 2
            else
                # 半基函数信息
                uniqueTri.isbd  =   true
                # 无负基函数部分
                faceii += 1
            end # if
        end #while

        #最后一个边是半基函数时进行处理
        faceii == 4*tetranum && begin
            # 挑出第swg2[faceii, 1]个基函数信息
            uniqueTri       =   uniqueTris[faceID]
            uniqueTri.isbd  =   true

        end # begin
    end
    # 构建 PWC 基函数信息
    pwcsInfo    =   [PWC{IntDtype, Precision.FT}() for _ in 1:3tetranum]
    # 按四面体编号排列，效果是相邻的四个三角形同属一个四面体并桉顺序排列
    face3opp1tri2unitri2   .=  sortslices(face3opp1tri2unitri2, dims = 1, by = x -> (x[5], x[6]), alg = ThreadsX.MergeSort)

    # 找出每个四面体对应的四个面
    for iTera in eachindex(tetrasInfo)
        # 该四面体
        tetra = tetrasInfo[iTera]
        # 四个三角形 id 
        tetra.faces       =     view(uniqueTris, unitri2[(4(iTera -1)+1):(4iTera), 1])
        # 四面体包含的基函数id
        resize!(tetra.inBfsID, 3)
        tetra.inBfsID       .=   (3(iTera -1)+1):(3iTera)
        # 写入 pwc 基函数信息
        for ii in 1:3
            pwcid   =   3(iTera - 1) + ii
            pwcsInfo[pwcid].bfID      =   pwcid
            pwcsInfo[pwcid].inGeo   =   iTera
            pwcsInfo[pwcid].center  =   tetra.center
        end
    end
    #tetraDataUnion(tetrasInfo, uniqueTris)
    return pwcsInfo

end





"""
计算构成六面体的所有三角形，并写入六面体包含的六个四边形的id，给分片常数 (PWC) 基函数赋值
"""
function setQuad4Hexas!(hexameshData::HexahedraMesh{IT, FT}, hexasInfo::Vector{HexahedraInfo{IT, FT, CT}}, ::Val{:PWC}) where {IT, FT, CT}
    # 六面体数
    hexnum      =   hexameshData.hexnum
    # 组成每个六面体的点的 id 数组
    hexasIDs    =   hexameshData.hexahedras
    # 组成每个面的局部 id 
    #* 注意不可改变下面的顺序因为跟计算基函数向量位置有关 *#
    quad6LocalIndices = quad6LocalIndices = [[2,3,7,6],[1,4,8,5],[4,3,7,8],[1,2,6,5],[5,6,7,8],[1,2,3,4]]

    #=  创建数组，用于后续保存所有面的点信息(4*2, 前四列排序
        寻找单独面，后四列保持四边形的点顺序)，然后是面所在的六面体id(1)，
        面所在 uniquad(独立四边形) 编号(1)、正负等信息(2)，共 12 列      =#
    face4hexd2bf2   =   zeros(IT, (hexnum*6, 12))
    # 利用内存视图提取其中的所有面的点
    face4           =   reshape(view(face4hexd2bf2, :, 1:4), (:, 24))
    face4quadids    =   reshape(view(face4hexd2bf2, :, 5:8), (:, 24))
    # 利用内存视图提取六面体信息的tri2和基函数信息的unitri2
    hex2            =   @view face4hexd2bf2[:, 9:10]
    uniquad2        =   @view face4hexd2bf2[:, 11:12]

    ## 创建六面体索引( id+localid )并根据排序更新
    # 此维度为六面体id
    hex2[:, 1]  .=   repeat(1:hexnum,  outer = 6 )
    # 此维度为该面在所在六面体的局部id
    hex2[:, 2]  .=   repeat(1:6,  inner = hexnum )
    
    # 将每个六面体的六个个面的点存放在 face4 数组中
    @threads for facei in 1:6
        viewofnodeids = view(hexasIDs, quad6LocalIndices[facei], :)'
        face4[:, facei:6:end] .=    viewofnodeids
        face4quadids[:, facei:6:end] .= viewofnodeids
    end
    # 重塑形状将每个边点组合放入同一维度
    face4        =  reshape(face4, (:, 4))
    face4quadids =  reshape(face4quadids, (:, 4))

    ##  对前四排序使得同一面的表示形式相同
    #   对前四个构成面的点排序
    @views sort!(face4hexd2bf2[:, 1:4], dims = 2)
    # 将每个六面体的四个角标点作为整体排序，次序按v1,v2,v3
    face4hexd2bf2   .=  sortslices(face4hexd2bf2, dims = 1, by = x -> (x[1], x[2], x[3], x[4]), alg = ThreadsX.MergeSort)

    ## for 循环计算PWC基函数id
    @inbounds let uniquadplusid   =   1
        uniquad2[1, 1]  =   1; uniquad2[1, 2] = 1
        # 初始值即1的位置要预计算一下
        uniquadplusid   +=  1

        # 从2开始循环，跟前一边不同为 unitri+，跟前一边不同为 unitri-
        for ii in 2:(hexnum*6)
            if  (face4[ii, 1]  != face4[ii-1, 1]) ||
                (face4[ii, 2]  != face4[ii-1, 2]) ||
                (face4[ii, 3]  != face4[ii-1, 3]) ||
                (face4[ii, 4]  != face4[ii-1, 4])
                uniquad2[ii, 1]      =   uniquadplusid
                uniquad2[ii, 2]      =   1
                uniquadplusid       +=   1
            elseif  (face4[ii, 1]  ==   face4[ii-1, 1]) &&
                (face4[ii, 2]  ==   face4[ii-1, 2]) && 
                (face4[ii, 3]  ==   face4[ii-1, 3]) && 
                (face4[ii, 4]  ==   face4[ii-1, 4]) 
                uniquad2[ii, 1]      =   uniquadplusid  - 1
                uniquad2[ii, 2]      =   -1
            end #if
        end #for
    end #begin

    # uniquadnum 函数数量( 同时也是独立四边形的数量 )
    uniquadnum   =   uniquad2[end, 1]
    # 单独三角形所在的 face4 的行 id
    uniquad2nodeRow    =   zeros(IT, uniquadnum)
    iunitri     =   1
    for inode in 1:(hexnum*6)
        (uniquad2[inode, 2] == 1) && begin
            uniquad2nodeRow[iunitri]     =   inode
            iunitri += 1
        end
    end
    # 预分配所有的构成六面体的四边形
    uniqueQuads  =   [Quads4Hexa{IT, FT}()  for _ in 1:uniquadnum ]

    # 设置 uniqueQuads 的三个点、边长、边法向量等信息
    @threads  for iQuad in 1:uniquadnum
        # 该四边形
        quad    =   uniqueQuads[iQuad]
        vertices    = quad.vertices

        # 该四边形对应的四个角点id
        quadvertIDs  =   view(face4quadids, uniquad2nodeRow[iQuad], :)
        # 设置角点
        for iv in 1:4
            vertices[:, iv]     .=   view(hexameshData.node, :, quadvertIDs[iv])
        end

        # 边指向向量
        edgev̂          =   quad.edgev̂
        @views for edgei in 1:4
            edgev̂[:, edgei]    =   vertices[:, edgei%4+1]  - vertices[:, edgei]
        end

        # 边长
        for edgei in 1:4
            @views quad.edgel[edgei]   =   norm(edgev̂[:, edgei])
        end

        @views facen̂    =  -cross(edgev̂[:, 1], edgev̂[:, 2])
        facen̂    ./=   norm(facen̂)
        # 计算边的指向单位向量、外法单位向量
        edgen̂                  =   quad.edgen̂
        for edgei in 1:4
            # 边的指向向量单位化
            @views edgev̂[:, edgei]  ./=   quad.edgel[edgei]
            # 最后计算边外法向量
            @views edgen̂[:, edgei]    =   cross(facen̂, edgev̂[:, edgei])
        end

    end

    # 写入公共面数是否为边界元信息
    @threads for faceii in 1:6hexnum
        # 面
        faceID   =   uniquad2[faceii, 1]
        ispOrm   =   uniquad2[faceii, 2]
        (ispOrm == -1) && begin
            uniqueQaud   =   uniqueQuads[faceID]
            uniqueQaud.isbd  =   false
        end
    end

    # 构造PWC
    pwcsInfo    =   [PWC{IntDtype, Precision.FT}() for _ in 1:3hexnum]
    # 按六面体编号排列，效果是相邻的四个四边形同属一个六面体并桉顺序排列
    face4hexd2bf2   .=  sortslices(face4hexd2bf2, dims = 1, by = x -> (x[9], x[10]), alg = ThreadsX.MergeSort)
    # 找出每个六面体对应的六个面
    @threads for iHex in eachindex(hexasInfo)
        # 该六面体
        hexa = hexasInfo[iHex]
        # 六个四边形 id 
        hexa.faces  .=     view(uniqueQuads, uniquad2[(6(iHex -1)+1):(6iHex), 1])
        # 六面体包含的基函数id
        resize!(hexa.inBfsID, 3)
        hexa.inBfsID       .=   (3(iHex -1)+1):(3iHex)
        # 写入 pwc 基函数信息
        for ii in 1:3
            pwcid   =   3(iHex - 1) + ii
            pwcsInfo[pwcid].bfID      =   pwcid
            pwcsInfo[pwcid].inGeo   =   iHex
            pwcsInfo[pwcid].center  =   hexa.center
        end
    end
    return pwcsInfo

end
