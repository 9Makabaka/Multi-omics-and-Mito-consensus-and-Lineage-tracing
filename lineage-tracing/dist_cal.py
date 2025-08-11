import numpy as np
from objects import DistObjects,RedeemR
import pandas as pd
import igraph as ig
from scipy.sparse import csr_matrix
import scipy.sparse as sp
from scipy.spatial.distance import pdist, squareform
from typing import List, Optional

def binary_dist(matrix: sp.csr_matrix, method: str) -> np.ndarray:
    """计算二元矩阵的距离指标（优化版）"""
    dense_matrix = matrix.toarray().astype(bool)
    n_cells = dense_matrix.shape[0]
    
    # 计算所有样本对的交集
    intersections = dense_matrix @ dense_matrix.T
    
    # 计算每个样本的元素和（修正：去掉.A1，改用NumPy的flatten()确保一维数组）
    sums = dense_matrix.sum(axis=1).flatten()  # 关键修改：用.flatten()替代.A1
    
    # 计算所有样本对的并集
    unions = sums[:, np.newaxis] + sums - intersections
    
    # 根据方法计算距离
    if method == "Jaccard":
        mask = unions == 0
        jaccard = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=np.float32), where=~mask)
        dist_matrix = 1 - jaccard
    elif method == "Dice":
        denominator = sums[:, np.newaxis] + sums
        mask = denominator == 0
        dice = np.divide(2 * intersections, denominator, out=np.zeros_like(intersections, dtype=np.float32), where=~mask)
        dist_matrix = 1 - dice
    elif method == "3WJaccard":
        mask = unions == 0
        jaccard3w = np.divide(intersections, unions, out=np.zeros_like(intersections, dtype=np.float32), where=~mask)
        dist_matrix = 1 - jaccard3w
    else:
        raise ValueError(f"不支持的距离方法: {method}")
    
    # 确保对角线为0
    np.fill_diagonal(dist_matrix, 0.0)
    return dist_matrix.astype(np.float32)


def quick_w_jaccard(matrix: sp.csr_matrix, weights: np.ndarray) -> np.ndarray:
    """计算加权Jaccard距离（优化版）"""
    dense_matrix = matrix.toarray().astype(bool)
    n_cells = dense_matrix.shape[0]
    
    # 应用权重
    weighted_matrix = dense_matrix * weights
    
    # 计算加权交集: (A * w) · B
    weighted_intersections = weighted_matrix @ dense_matrix.T
    
    # 计算每个样本的加权和
    weighted_sums = weighted_matrix.sum(axis=1)
    
    # 计算加权并集: sum(A*w) + sum(B*w) - sum(min(A,B)*w)
    weighted_unions = weighted_sums[:, np.newaxis] + weighted_sums - weighted_intersections
    
    # 计算加权Jaccard距离
    mask = weighted_unions == 0
    w_jaccard = np.divide(weighted_intersections, weighted_unions, 
                         out=np.zeros_like(weighted_intersections, dtype=np.float32), 
                         where=~mask)
    dist_matrix = 1 - w_jaccard
    
    # 确保对角线为0
    np.fill_diagonal(dist_matrix, 0.0)
    return dist_matrix.astype(np.float32)


def quick_w_cosine(matrix: sp.csr_matrix, weights: np.ndarray) -> np.ndarray:
    """计算加权余弦距离（优化版）"""
    dense_matrix = matrix.toarray()
    
    # 应用权重
    weighted_matrix = dense_matrix * weights
    
    # 计算点积
    dot_products = weighted_matrix @ weighted_matrix.T
    
    # 计算每个向量的范数
    norms = np.linalg.norm(weighted_matrix, axis=1, keepdims=True)
    norm_products = norms @ norms.T
    
    # 计算余弦相似度
    mask = norm_products == 0
    cos_sim = np.divide(dot_products, norm_products, 
                       out=np.zeros_like(dot_products, dtype=np.float32), 
                       where=~mask)
    
    # 余弦距离 = 1 - 余弦相似度
    dist_matrix = 1 - cos_sim
    
    # 确保对角线为0
    np.fill_diagonal(dist_matrix, 0.0)
    return dist_matrix.astype(np.float32)


def add_dist(
    redeemr_obj: RedeemR,  
    jaccard: bool = True,
    dice: bool = True,
    jaccard3w: bool = True,
    w_jaccard: bool = True,
    w_cosine: bool = True,
    weight_df: Optional[pd.DataFrame] = None,
    nn: int = 1,
    lsi_dist: bool = True,
    dim: List[int] = list(range(1, 50))
) -> RedeemR:
    """计算并添加多种距离指标（优化版）"""
    # 初始化距离对象（假设DistObjects已定义）
    dist_objects = DistObjects()
    
    # 检查二值化矩阵是否存在
    if redeemr_obj.cts_mtx_bi is None:
        raise ValueError("请先运行make_matrix生成二值化矩阵(cts_mtx_bi)")
    
    binary_matrix = redeemr_obj.cts_mtx_bi["matrix"]
    variants = redeemr_obj.cts_mtx_bi["columns"]
    
    # 处理权重数据（与原代码相同）
    weight = None
    if weight_df is not None:
        if "Variant" not in weight_df.columns or "CellPCT" not in weight_df.columns:
            raise ValueError("weight_df必须包含'Variant'和'CellPCT'列")
        
        weight_df = weight_df.copy()
        weight_df["weight"] = 1 - weight_df["CellPCT"]
        
        variant_weights = pd.DataFrame({"Variant": variants})
        print("variant_weights 'Variant' dtype:", variant_weights["Variant"].dtype)
        print("weight_df 'Variant' dtype:", weight_df["Variant"].dtype)

        # 两列统一转换为字符串类型（object）
        variant_weights["Variant"] = variant_weights["Variant"].astype(str)
        weight_df["Variant"] = weight_df["Variant"].astype(str)
        
        variant_weights = variant_weights.merge(weight_df[["Variant", "weight"]], 
                                              on="Variant", how="left")
        
        na_count = variant_weights["weight"].isna().sum()
        if na_count > 0:
            variant_weights["weight"].fillna(nn, inplace=True)
            print(f"注意：{na_count}个变异在权重数据中未找到，使用默认权重{nn}")
        
        if len(variant_weights) != len(variants):
            raise ValueError(f"权重长度({len(variant_weights)})与变异数量({len(variants)})不匹配")
        
        weight = variant_weights["weight"].values
        print("权重向量与细胞-变异矩阵匹配，继续计算...")
    
    # 计算各种距离（使用优化后的函数）
    if jaccard:
        dist_objects.jaccard = binary_dist(binary_matrix, "Jaccard")
        print("已添加jaccard距离")
    
    if dice:
        dist_objects.dice = binary_dist(binary_matrix, "Dice")
        print("已添加dice距离")
    
    if jaccard3w:
        dist_objects.jaccard3W = binary_dist(binary_matrix, "3WJaccard")
        print("已添加3wjaccard距离")
    
    if w_jaccard:
        if weight is None:
            raise ValueError("计算加权jaccard距离需要提供weight_df参数")
        dist_objects.w_jaccard = quick_w_jaccard(binary_matrix, weight)
        print("已添加weighted jaccard距离")
    
    if w_cosine:
        if weight is None:
            raise ValueError("计算加权cosine距离需要提供weight_df参数")
        dist_objects.w_cosine = quick_w_cosine(binary_matrix, weight)
        print("已添加weighted cosine距离")
    
    if lsi_dist:
        if redeemr_obj.seurat is None or "X_pca_subset" not in redeemr_obj.seurat.obsm:
            raise ValueError("请先运行seurat_lsi_clustering生成LSI结果")
        
        lsi_embeddings = redeemr_obj.seurat.obsm["X_pca_subset"]
        # 验证维度数量是否正确（应等于len(dim)=49）
        if lsi_embeddings.shape[1] != len(dim):
            raise ValueError(f"LSI嵌入维度不匹配，实际{lsi_embeddings.shape[1]}，预期{len(dim)}")
    
        dist_objects.lsi_dist = squareform(pdist(lsi_embeddings, metric="euclidean"))
        print("已添加LSI距离")
    
    redeemr_obj.dist_objects = dist_objects
    return redeemr_obj

def from_dist2graph(d, k_param=30, return_igraph=True):
    """
    从距离矩阵生成互近邻(MNN)网络，对应R中的FromDist2Graph函数
    
    参数:
        d: 距离矩阵（可以是numpy数组或类似矩阵的结构），行和列均为细胞
        k_param: 每个细胞需要寻找的最近邻数量，默认30
        return_igraph: 是否返回igraph对象，默认True（返回igraph网络），False则返回邻接矩阵
    
    返回:
        igraph.Graph或scipy.sparse.csr_matrix: 互近邻网络（igraph对象）或邻接矩阵
    """
    # 1. 确保距离矩阵为numpy数组格式
    if not isinstance(d, np.ndarray):
        try:
            d = np.array(d)  # 尝试转换为numpy数组
        except:
            raise TypeError("距离矩阵必须可转换为numpy数组")
    
    # 检查矩阵是否为方阵
    if d.shape[0] != d.shape[1]:
        raise ValueError("距离矩阵必须是方阵（行和列数相等）")
    
    n_cells = d.shape[0]  # 细胞数量
    cell_names = np.arange(n_cells)  # 默认为索引，可根据实际情况替换为细胞ID
    
    # 2. 计算每个细胞的K最近邻（KNN）索引
    # 对每行距离排序，取前k_param个索引（排除自身，若存在）
    knn_indices = np.zeros((n_cells, k_param), dtype=int)
    for i in range(n_cells):
        # 对距离排序并获取索引（从小到大）
        sorted_indices = np.argsort(d[i, :])
        # 排除自身（如果第一个是自己）
        if sorted_indices[0] == i:
            knn_indices[i, :] = sorted_indices[1:k_param+1]  # 取第2到第k_param+1个
        else:
            knn_indices[i, :] = sorted_indices[:k_param]  # 取前k_param个
    
    # 3. 筛选互近邻（MNN）：仅保留双向近邻关系
    edges = []
    for i in range(n_cells):
        # 当前细胞i的K个近邻
        neighbors = knn_indices[i, :]
        # 检查这些近邻是否也将i视为近邻
        mutual_mask = [i in knn_indices[neigh, :] for neigh in neighbors]
        # 获取互近邻索引
        mutual_neighbors = neighbors[mutual_mask]
        
        # 如果没有互近邻，至少保留与最近的1个细胞的连接
        if len(mutual_neighbors) == 0:
            mutual_neighbors = [neighbors[0]]  # 取最近的一个
        
        # 添加边（确保i < j，避免重复边）
        for j in mutual_neighbors:
            if i < j:  # 只添加一次，避免双向重复
                edges.append((i, j))
    
    # 4. 构建邻接矩阵（稀疏矩阵）
    # 提取边的行和列索引
    rows, cols = zip(*edges) if edges else ([], [])
    # 构建稀疏邻接矩阵（值为1表示存在连接）
    adj_matrix = csr_matrix(
        (np.ones(len(rows), dtype=int), (rows, cols)),
        shape=(n_cells, n_cells)
    )
    # 确保矩阵对称（无向图）
    adj_matrix = adj_matrix + adj_matrix.T
    adj_matrix[adj_matrix > 1] = 1  # 避免重复边导致的值大于1
    
    # 5. 设置细胞名称（行和列名）
    adj_matrix = adj_matrix.tocsr()
    adj_matrix.cell_names = cell_names  # 附加细胞名称属性
    
    # 6. 返回结果
    if return_igraph:
        # 转换为igraph对象
        g = ig.Graph.Adjacency(adj_matrix.toarray().tolist(), mode="undirected")
        # 设置节点名称
        g.vs["name"] = [str(name) for name in cell_names]
        return g
    else:
        return adj_matrix
    