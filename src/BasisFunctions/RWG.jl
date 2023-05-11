# 定义RWG基函数数据结构
"""
    RWG{IT<:Integer , FT<:AbstractFloat} <: LinearBasisFunction

RWG基函数复合类型：
```
isbd        ::Bool              是否为边界元即半基函数，布尔类型
bfID        ::IT                基函数编号，整形
edgel       ::FT                基函数边长，浮点型
inGeo       ::MVector{2, IT}    基函数所在两个三角形编号（半基函数为1个，赋值成一样的两个），长度为2的向量数组
inGeoID     ::MVector{2, IT}    基函数在两个上面三角形中的局部编号（半基函数为1个，赋值成一样的两个），取值1:3，长度为2的向量数组
center      ::MVec3D{FT}        基函数中心，用于八叉树分组
```
"""
mutable struct RWG{IT<:Integer , FT<:AbstractFloat} <: LinearBasisFunction
    isbd        ::Bool
    bfID        ::IT
    edgel       ::FT
    inGeo       ::MVector{2, IT}
    inGeoID     ::MVector{2, IT}
    center      ::MVec3D{FT}
end

"""
    RWG{IT, FT}()where {IT <: Integer, FT<:AbstractFloat}

RWG 的默认构造函数，默认非面的边界，所有元素置零。
"""
function RWG{IT, FT}()where {IT <: Integer, FT<:AbstractFloat}
    isbd        =    false
    bfID        =    zero(IT)  
    edgel       =    zero(FT)
    inGeo       =    zero(MVector{2, IT})
    inGeoID     =    zero(MVector{2, IT})
    center      =    zero(MVec3D{FT})
    return RWG{IT, FT}(isbd, bfID,  edgel,  inGeo, inGeoID, center)
end
RWG()   =    RWG{IntDtype, Precision.FT}()

"""
    rwgbfConstructerTrianglesInfoModifiers!(trianglemeshData::TriangleMesh, trianglesInfo::Vector{TriangleInfo{IT, FT}}) where {IT<:Integer, FT<:AbstractFloat}

此函数采用排序算法，将每个边的两个点、在三角形中的对点、三角形、基函数等属性信息放在一个大数组中，
通过对不同的属性排序（如按边所在点随数组排序即可将边相同的点放在一起），即可得到基函数分组信息，以此可构造RWG基函数。
函数完成以下功能：
1. 构造基函数类型实例数组（类似结构化数组）rwgsInfo记录基函数相关信息，并作为返回值；
2. 写入三角形类型实例数组trianglesInfo中关于基函数的信息。
"""
function rwgbfConstructerTrianglesInfoModifiers!(trianglemeshData::TriangleMesh, trianglesInfo::Vector{TriangleInfo{IT, FT}}) where {IT<:Integer, FT<:AbstractFloat}
    # 三角形数
    trinum  =   trianglemeshData.trinum
    
    # 创建数组，用于后续保存所有边的点和边的对点信息，然后是边所在的三角形id，在三角形中是第几个边，边所在RWG编号、正负等信息，共6列
    edge2opp1tri2rwg2   =   zeros(IT, (trinum*3, 7))
    # 利用内存视图提取其中的所有边的点和边的对点信息信息
    edge2opp1           =   reshape(view(edge2opp1tri2rwg2, :, 1:3), (:, 9))
    # 利用内存视图提取三角形信息的tri2和基函数信息的rwg2
    tri2                =   @view edge2opp1tri2rwg2[:, 4:5]
    rwg2                =   @view edge2opp1tri2rwg2[:, 6:7]
    
    ## 创建三角形索引并根据排序更新
    # 此维度为
    tri2[:, 1] =   repeat(1:trinum,     outer=3)
    tri2[:, 2] =   repeat(1:3,  inner = trinum )
    
    # 将每个三角形的三条边的点及其对点存放在edge2opp1数组中
    #* 注意不可改变下面的顺序因为跟RWG基函数对应的自由点在三角形中顺序有关 *#
    edge2opp1[:, 1] =   edge2opp1[:, 6] =   edge2opp1[:, 8] =   trianglemeshData.triangles[2, :]
    edge2opp1[:, 2] =   edge2opp1[:, 4] =   edge2opp1[:, 9] =   trianglemeshData.triangles[3, :]
    edge2opp1[:, 3] =   edge2opp1[:, 5] =   edge2opp1[:, 7] =   trianglemeshData.triangles[1, :]
    # 重塑形状将每个边点组合放入同一维度
    edge2opp1        =   reshape(edge2opp1, (:, 3))

    ##  对前两列排序使得同一边的表示形式相同
    #   对前两个构成边的点排序
    @views sort!(edge2opp1tri2rwg2[:, 1:2], dims = 2, alg = ThreadsX.MergeSort)

    # 将每个三角形的三个角标点作为整体排序，次序按v1,v2,v3
    edge2opp1tri2rwg2[:]   =  sortslices(edge2opp1tri2rwg2, dims = 1, by = x -> (x[1], x[2]), alg = ThreadsX.MergeSort)

    ## for 循环计算RWG基函数id
    rwgplusid = @inbounds let rwgplusid   =   1
        rwg2[1, 1]  =   1; rwg2[1, 2] = 1
        # 初始值即1的位置要预计算一下
        rwgplusid   +=  1

        # 从2开始循环，跟前一边不同为RWG+，跟前一边不同为RWG-
        for ii in 2:(trinum*3)
            if  (edge2opp1[ii, 1]  != edge2opp1[ii-1, 1]) ||
                (edge2opp1[ii, 2]  != edge2opp1[ii-1, 2])
                rwg2[ii, 1]     =   rwgplusid
                rwg2[ii, 2]     =   1
                rwgplusid   +=  1
            elseif  (edge2opp1[ii, 1]  ==   edge2opp1[ii-1, 1]) &&
                (edge2opp1[ii, 2]  ==   edge2opp1[ii-1, 2])
                rwg2[ii, 1]        =   rwgplusid  - 1
                rwg2[ii, 2]        =   -1
            end #if
        end #for

        rwgplusid
    end #begin

    # rwg函数数量
    rwgnum      =   rwgplusid - 1
    # 保存RWG基函数信息的数组
    rwgsInfo    =   [RWG() for _ in 1:rwgnum]

    # 写入基函数信息
    @views let edgeii = 1, bfID::IT
        while edgeii < 3*trinum
            bfID   =   rwg2[edgeii, 1]
            # 挑出第2[edgeii,1]个基函数信息
            rwgInfo =   rwgsInfo[bfID]
            if rwg2[edgeii, 1] ==  rwg2[edgeii + 1, 1]
                # 全基函数信息
                rwgInfo.isbd        =   false
                rwgInfo.bfID       =   bfID
                rwgInfo.inGeo      .=   tri2[edgeii:edgeii+1, 1]
                rwgInfo.inGeoID    .=   tri2[edgeii:edgeii+1, 2]
                # 设置边所在的点、计算边长
                edgenode            =   trianglemeshData.node[:,edge2opp1[edgeii,1:2]]
                rwgInfo.edgel       =   norm(edgenode[:, 1] - edgenode[:, 2])
                rwgInfo.center[:]   =   mean(edgenode, dims = 2)
                # 跳过负基函数部分
                edgeii += 2
            else
                # 半基函数信息
                rwgInfo.isbd    =   true
                rwgInfo.bfID   =   bfID
                rwgInfo.inGeo[1]=   tri2[edgeii, 1]
                rwgInfo.inGeoID[1]  =   tri2[edgeii, 2]
                # 设置边所在的点、计算边长
                edgenode            =   trianglemeshData.node[:,edge2opp1[edgeii,1:2]]
                rwgInfo.edgel       =   norm(edgenode[:, 1] - edgenode[:, 2])
                rwgInfo.center[:]   =   mean(edgenode, dims = 2)
                # 无负基函数部分
                edgeii += 1
            end # if
        end #while

        #最后一个边是半基函数时进行处理
        edgeii == 3*trinum && begin
            # 挑出第rwg2[edgeii, 1]个基函数信息
            rwgInfo         =   rwgsInfo[rwgnum]
            rwgInfo.isbd    =   true
            rwgInfo.bfID   =   bfID
            rwgInfo.inGeo   =   tri2[edgeii, 1]
            rwgInfo.inGeoID =   tri2[edgeii, 2]
            edgenode        =   trianglemeshData.node[:,edge2opp1[edgeii,1:2]]
            rwgInfo.edgel   =   norm(edgenode[:, 1] - edgenode[:, 2])
            rwgInfo.center[:]   =   mean(edgenode, dims = 2)
        end # begin

    end #let

    ## 设置三角形所包含的基函数信息
    #  将edge2opp1tri2rwg2按三角形序号重新排序，以将每个三角形包含的基函数信息放在一起
    edge2opp1tri2rwg2      .=  sortslices(edge2opp1tri2rwg2, dims = 1, by = x -> (x[4], x[5]), alg = ThreadsX.MergeSort)
    
    #   此时bfid\正负信息每连续3个同属一个三角形，可将三角形包含的基函数信息写入三角形信息中
    @views rwgsID  =   reshape(rwg2[:, 1], (3,:))
    @views rwgspm  =   reshape(rwg2[:, 2], (3,:))
    for triidx in 1:length(trianglesInfo)
        triInfo =   trianglesInfo[triidx]
        triInfo.inBfsID     .=   view(rwgsID, :, triidx)

        # 此处将rwg基函数的正负直接封装进边长的正负
        for edgeii in 1:3
            (rwgspm[edgeii,triidx] == -1) &&  (triInfo.edgel[edgeii] *= -1)
        end
    end

    return rwgsInfo

end #function

"""
    rwgbfnohalfConstructerTrianglesInfoModifiers!(trianglemeshData::TriangleMesh, trianglesInfo::Vector{TriangleInfo{IT, FT}}) where {IT<:Integer, FT<:AbstractFloat}

此函数与[`rwgbfConstructerTrianglesInfoModifiers!`](@ref)基本一致，不同在于不生成半基函数。
"""
function rwgbfnohalfConstructerTrianglesInfoModifiers!(trianglemeshData::TriangleMesh, trianglesInfo::Vector{TriangleInfo{IT, FT}}) where {IT<:Integer, FT<:AbstractFloat}
    # 三角形数
    trinum  =   trianglemeshData.trinum
    
    # 创建数组，用于后续保存所有边的点和边的对点信息，然后是边所在的三角形id，在三角形中是第几个边，边所在RWG编号、正负等信息，共6列
    edge2opp1tri2rwg2   =   zeros(IT, (trinum*3, 7))
    # 利用内存视图提取其中的所有边的点和边的对点信息信息
    edge2opp1           =   reshape(view(edge2opp1tri2rwg2, :, 1:3), (:, 9))
    # 利用内存视图提取三角形信息的tri2和基函数信息的rwg2
    tri2                =   @view edge2opp1tri2rwg2[:, 4:5]
    rwg2                =   @view edge2opp1tri2rwg2[:, 6:7]
    
    ## 创建三角形索引并根据排序更新
    # 此维度为
    tri2[:, 1] =   repeat(1:trinum,     outer=3)
    tri2[:, 2] =   repeat(1:3,  inner = trinum )
    
    # 将每个三角形的三条边的点及其对点存放在edge2opp1数组中
    #* 注意不可改变下面的顺序因为跟RWG基函数对应的自由点在三角形中顺序有关 *#
    edge2opp1[:, 1] =   edge2opp1[:, 6] =   edge2opp1[:, 8] =   trianglemeshData.triangles[2, :]
    edge2opp1[:, 2] =   edge2opp1[:, 4] =   edge2opp1[:, 9] =   trianglemeshData.triangles[3, :]
    edge2opp1[:, 3] =   edge2opp1[:, 5] =   edge2opp1[:, 7] =   trianglemeshData.triangles[1, :]
    # 重塑形状将每个边点组合放入同一维度
    edge2opp1        =   reshape(edge2opp1, (:, 3))

    ##  对前两列排序使得同一边的表示形式相同
    #   对前两个构成边的点排序
    @views sort!(edge2opp1tri2rwg2[:, 1:2], dims = 2)

    # 将每个三角形的三个角标点作为整体排序，次序按v1,v2,v3
    edge2opp1tri2rwg2[:]   =  sortslices(edge2opp1tri2rwg2, dims = 1, by = x -> (x[1], x[2]), alg = ThreadsX.MergeSort)

    ## for 循环计算RWG基函数id
    rwgplusid = @inbounds let rwgplusid   =   1
        # 从1开始循环，跟后一边相同同为RWG+，跟前一边不同为RWG-
        for ii in 1:(trinum*3 - 1)
            if  (edge2opp1[ii, 1]  == edge2opp1[ii+1, 1]) &&
                (edge2opp1[ii, 2]  == edge2opp1[ii+1, 2])
                # 正基函数
                rwg2[ii, 1]     =   rwgplusid
                rwg2[ii, 2]     =   1
                # 负基函数
                rwg2[ii+1, 1]   =   rwgplusid
                rwg2[ii+1, 2]   =  -1

                rwgplusid   +=  1
            end #if
        end #for
        rwgplusid
    end #begin

    # rwg函数数量
    rwgnum      =   rwgplusid - 1
    # 保存RWG基函数信息的数组
    rwgsInfo    =   [RWG() for _ in 1:rwgnum]

    # 写入基函数信息
    @views let edgeii = 1, bfID::IT
        for edgeii in 1:(trinum*3 - 1)
            bfID   =   rwg2[edgeii, 1]
            # id是0则为半基函数，此处不考虑，跳过
            bfID   == 0  && continue
            # 跳过－基函数
            rwg2[edgeii, 2] ==  -1 && continue

            # 挑出第2[edgeii,1]个基函数信息
            rwgInfo =   rwgsInfo[bfID]
            if rwg2[edgeii, 1] ==  rwg2[edgeii + 1, 1]
                # 全基函数信息
                rwgInfo.isbd        =   false
                rwgInfo.bfID       =   bfID
                rwgInfo.inGeo      .=   tri2[edgeii:edgeii+1, 1]
                rwgInfo.inGeoID    .=   tri2[edgeii:edgeii+1, 2]
                # 设置边所在的点、计算边长
                edgenode            =   trianglemeshData.node[:,edge2opp1[edgeii,1:2]]
                rwgInfo.edgel       =   norm(edgenode[:, 1] - edgenode[:, 2])
                rwgInfo.center[:]   =   mean(edgenode, dims = 2)
            end # if
        end #while
    end #let

    ## 设置三角形所包含的基函数信息
    #  将edge2opp1tri2rwg2按三角形序号重新排序，以将每个三角形包含的基函数信息放在一起
    edge2opp1tri2rwg2[:]   =  sortslices(edge2opp1tri2rwg2, dims = 1, by = x -> (x[4], x[5]), alg = ThreadsX.MergeSort)
    
    #   此时bfid\正负信息每连续3个同属一个三角形，可将三角形包含的基函数信息写入三角形信息中
    @views rwgsID  =   reshape(rwg2[:, 1], (3,:))
    @views rwgspm  =   reshape(rwg2[:, 2], (3,:))
    @views for triidx in 1:length(trianglesInfo)
        triInfo =   trianglesInfo[triidx]
        triInfo.inBfsID     =   rwgsID[:,triidx]

        # 此处将rwg基函数的正负直接封装进边长的正负
        for edgeii in 1:3
            @views (rwgspm[edgeii,triidx] == -1) &&  (triInfo.edgel[edgeii] *= -1)
        end
    end

    return rwgsInfo

end #function

"""
    getTriangleInfo(trianglemeshData::TriangleMesh{IT, FT}) where{IT, FT}

根据网格信息 `trianglemeshData` 生成三角形信息 `trianglesInfo` 、RWG基函数信息 `rwgsInfo`。
"""
function getTriangleInfo(trianglemeshData::TriangleMesh{IT, FT}) where{IT, FT}
    # 分配三角形信息类型数组
    trianglesInfo       =   [TriangleInfo{IT, FT}()  for _ in 1:trianglemeshData.trinum ]
    # 写入对应的三角形id和点坐标数据
    setTricoor!(trianglesInfo, trianglemeshData)
    # 写入三角形边长、边外法向量、面法向量、面积
    setTriParam!.(trianglesInfo)
    # 计算基函数信息并修改三角形信息中基函数相关项
    rwgsInfo    =   rwgbfnohalfConstructerTrianglesInfoModifiers!(trianglemeshData, trianglesInfo)
    return trianglesInfo, rwgsInfo
end
