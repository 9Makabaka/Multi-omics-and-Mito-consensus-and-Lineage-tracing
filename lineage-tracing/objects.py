from typing import Dict, List, Optional
import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
from scipy.spatial.distance import squareform
import scanpy as sc
from sklearn.feature_extraction.text import TfidfTransformer
from typing import Optional, List, Dict, Union
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.BaseTree import Tree


class DistObjects:
    """存储各种距离计算结果的容器类，对应R中的DistObjects"""
    def __init__(self):
        self.jaccard: Optional[np.ndarray] = None          # Jaccard距离
        self.dice: Optional[np.ndarray] = None            # Dice距离
        self.jaccard3W: Optional[np.ndarray] = None       # 3W Jaccard距离
        self.w_jaccard: Optional[np.ndarray] = None       # 加权Jaccard距离
        self.w_cosine: Optional[np.ndarray] = None        # 加权余弦距离
        self.lsi_dist: Optional[np.ndarray] = None        # LSI距离
        
class TREE:
    """模拟R中的TREE类，存储系统发育树相关数据"""
    def __init__(self):
        self.phylo: Optional[Tree] = None  # Bio.Phylo的Tree对象
        self.treedata: Optional[Dict] = None  # 树的结构化数据
        self.records: Optional[str] = None  # 记录构建参数

class RedeemR:
    """模拟R中的redeemR类，存储mtDNA突变分析的结构化数据"""
    def __init__(self):
        # 过滤后的基因型数据（每个细胞的突变记录）
        self.gt_summary_filtered: Optional[pd.DataFrame] = None
        # 合格细胞的元数据（细胞ID及覆盖度等）
        self.cell_meta: Optional[pd.DataFrame] = None
        # 过滤后的变异特征（每个变异的汇总统计）
        self.v_filtered: Optional[pd.DataFrame] = None
        # 纯质性变异列表
        self.homo_variants: Optional[List[str]] = None
        # 唯一变异列表
        self.unique_v: Optional[List[str]] = None
        # 深度信息摘要（细胞覆盖度等）
        self.depth_summary: Optional[Dict] = None
        # 分析参数（阈值、切割值等）
        self.para: Optional[Dict] = None
        # 附加信息（过滤记录、数据路径等）
        self.attr: Optional[Dict] = None
        # 新增：细胞-变异计数矩阵（稀疏矩阵）
        self.cts_mtx: Optional[Dict] = None  # {matrix: csr_matrix, rows: cell_ids, columns: variants}
        # 新增：二值化矩阵（稀疏矩阵）
        self.cts_mtx_bi: Optional[Dict] = None  # 同上
        # 新增：存储单细胞分析结果（替代Seurat对象）
        self.seurat: Optional[ad.AnnData] = None
        # 新增：距离计算结果
        self.dist_objects: Optional[DistObjects] = None
        # 新增：系统发育树
        self.tree: Optional[TREE] = None

    def __repr__(self) -> str:
        return f"RedeemR object with {len(self.unique_v) if self.unique_v else 0} unique variants"
    def seurat_lsi_clustering(self, binary: bool = True, res: float = 0.6, 
                             lsi_dim: List[int] = list(range(1, 50)),  # 对应R中的2:50（Python索引从0开始）
                             rm_variants: List[str] = ["Variants310TC", "Variants3109TC", "Variants5764CT"]):
        """
        使用线粒体变异进行基于LSI的聚类分析，对应R中的SeuratLSIClustering
        
        参数:
            binary: 是否使用二值化矩阵
            res: 聚类分辨率（Leiden算法）
            lsi_dim: 用于聚类的LSI维度（索引从0开始）
            rm_variants: 需要移除的变异列表
        
        返回:
            自身对象（更新后的RedeemR实例）
        """
        # 1. 选择使用的矩阵（二值化或计数矩阵）
        if binary:
            if self.cts_mtx_bi is None:
                raise ValueError("二值化矩阵(cts_mtx_bi)不存在，请先调用make_matrix生成")
            matrix_data = self.cts_mtx_bi["matrix"].copy()
            variants = self.cts_mtx_bi["columns"].copy()
            cells = self.cts_mtx_bi["rows"].copy()
        else:
            if self.cts_mtx is None:
                raise ValueError("计数矩阵(cts_mtx)不存在，请先调用make_matrix生成")
            matrix_data = self.cts_mtx["matrix"].copy()
            variants = self.cts_mtx["columns"].copy()
            cells = self.cts_mtx["rows"].copy()
        
        # 2. 移除指定变异
        keep_var_idx = [i for i, var in enumerate(variants) if var not in rm_variants]
        if not keep_var_idx:
            raise ValueError("所有变异均被移除，请检查rm_variants参数")
        matrix_data = matrix_data[:, keep_var_idx]
        variants = [variants[i] for i in keep_var_idx]
        
        # 3. 过滤无突变的细胞（行和为0的细胞）
        cell_sums = matrix_data.sum(axis=1).A.flatten() if sp.sparse.issparse(matrix_data) else matrix_data.sum(axis=1)
        keep_cell_idx = cell_sums > 0
        if not np.any(keep_cell_idx):
            raise ValueError("所有细胞均无突变，无法进行聚类分析")
        matrix_data = matrix_data[keep_cell_idx, :]
        cells = [cells[i] for i, keep in enumerate(keep_cell_idx) if keep]
        
        # 4. 创建AnnData对象（细胞为行，变异为列）
        adata = ad.AnnData(
            X=matrix_data,  # 细胞×变异矩阵（行=细胞，列=变异）
            obs=pd.DataFrame(index=cells),  # 细胞元数据（行数=细胞数）
            var=pd.DataFrame(index=variants)  # 变异元数据（列数=变异数）
        )
        # 设置所有变异为高度可变特征（对应原R代码的VariableFeatures设置）
        adata.var['highly_variable'] = True
    
        # 5. TF-IDF转换（对应R中的RunTFIDF）
        tfidf = TfidfTransformer(norm="l2")
        adata.X = tfidf.fit_transform(adata.X)  # 对细胞×变异矩阵进行转换
    
        # # 6. 筛选至少在2个细胞中出现的变异（对应R中的FindTopFeatures）
        # sc.pp.filter_genes(adata, min_cells=2)  # 过滤变异（列）
    
        # 7. 执行PCA（替代LSI，原R代码用RunSVD）
        sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
        # lsi_dim是需要使用的维度索引（如[1,2,...,49]）
        adata.obsm['X_pca_subset'] = adata.obsm['X_pca'][:, lsi_dim]  # 提取子空间
    
        # 8. 非线性降维（UMAP和t-SNE）
        # 使用提取的子空间计算邻居，通过use_rep指定
        sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca_subset')  # 替代dims参数
        sc.tl.umap(adata)
        sc.tl.tsne(adata, use_rep='X_pca_subset')  # 基于子空间降维
    
        # 9. 聚类分析（Leiden算法）
        sc.tl.leiden(adata, resolution=res)  # 结果存储在obs['leiden']
    
        # 10. 存储结果
        self.seurat = adata
        print(f"LSI聚类完成：{len(adata.obs['leiden'].unique())}个聚类")
    
        return self
    
    def make_tree(self, d: str, algorithm: str, only_return_tree: bool = False) -> Union['RedeemR', TREE]:
        """
        生成系统发育树
        
        参数:
            d: 距离类型，可选"jaccard"、"Dice"、"jaccard3W"、"w_jaccard"、"w_cosine"、"LSIdist"
            algorithm: 构建树的算法，可选"nj"（邻接法）或"upgma"
            only_return_tree: 是否只返回树对象，而非完整的RedeemR对象
        
        返回:
            更新后的RedeemR对象或TREE对象
        """
        # 验证距离类型是否有效
        valid_dist_types = ["jaccard", "Dice", "jaccard3W", "w_jaccard", "w_cosine", "LSIdist"]
        if d not in valid_dist_types:
            raise ValueError(f"无效的距离类型: {d}，必须是{valid_dist_types}之一")
        
        # 验证算法是否有效
        if algorithm not in ["nj", "upgma"]:
            raise ValueError(f"无效的算法: {algorithm}，必须是'nj'或'upgma'")
        
        # 验证距离矩阵是否存在
        if not self.dist_objects:
            raise ValueError("dist_objects未初始化，请先计算距离矩阵")
        
        # 获取指定的距离矩阵
        dist_matrix = getattr(self.dist_objects, d, None)
        if dist_matrix is None:
            raise ValueError(f"距离矩阵 {d} 不存在于dist_objects中")
        
        # 准备距离矩阵（转换为Bio.Phylo所需的格式）
        # 假设距离矩阵是方阵，且细胞ID在cts_mtx的rows中
        cell_ids = self.cts_mtx['rows'] if self.cts_mtx else [f"cell_{i}" for i in range(dist_matrix.shape[0])]
        
        # 将方阵转换为三角矩阵（Bio.Phylo要求的格式）
        try:
            # 1. 验证距离矩阵是对称的数值方阵（同之前）
            if not np.allclose(dist_matrix, dist_matrix.T):
                raise ValueError("距离矩阵必须是对称矩阵")
            if dist_matrix.dtype.kind not in 'iufc':
                raise TypeError("距离矩阵必须包含数值类型元素")
            n = dist_matrix.shape[0]
            if dist_matrix.shape[1] != n:
                raise ValueError(f"距离矩阵必须是方阵，实际形状：{dist_matrix.shape}")
            if len(cell_ids) != n:
                cell_ids = cell_ids[:n]
                print(f"修正 cell_ids 长度为 {n}（与距离矩阵匹配）")

            # 2. 生成符合 "每行长度1,2,...,n" 的下三角矩阵（包含对角线0）
            dist_list = []
            for i in range(n):
                # 提取第i行的前i个元素（到样本0~i-1的距离）
                row_without_diag = dist_matrix[i, :i].tolist()
                # 追加对角线元素（自身到自身的距离0）
                row_with_diag = row_without_diag + [0.0]  # 关键：添加对角线0
                # 校验每行长度是否为 i+1（符合 range(1, n+1)）
                if len(row_with_diag) != i + 1:
                    raise ValueError(f"第{i}行格式错误：应有{i+1}个元素，实际{len(row_with_diag)}个")
                dist_list.append(row_with_diag)

            # 3. 校验整体格式（同之前）
            if len(dist_list) != len(cell_ids):
                raise ValueError(f"矩阵行数（{len(dist_list)}）与样本数（{len(cell_ids)}）不匹配")
            if [len(row) for row in dist_list] != list(range(1, len(cell_ids) + 1)):
                raise ValueError(f"矩阵格式不符合要求，每行长度应为{list(range(1, len(cell_ids)+1))}")

            # 4. 创建 DistanceMatrix 对象
            dm = DistanceMatrix(names=cell_ids, matrix=dist_list)
        except Exception as e:
            raise RuntimeError(f"距离矩阵转换失败: {str(e)}")
        
        # 构建系统发育树
        try:
            constructor = DistanceTreeConstructor()
            if algorithm == "nj":
                phylo_tree = constructor.nj(dm)  # 使用nj方法
            else:  # upgma
                phylo_tree = constructor.upgma(dm)  # 使用upgma方法
        except Exception as e:
            raise RuntimeError(f"构建系统发育树失败: {str(e)}")
        
        # 创建TREE对象
        tree_obj = TREE()
        tree_obj.phylo = phylo_tree
        tree_obj.records = f"{d}-{algorithm}"
        
        # 转换为treedata格式（简化版，存储节点和分支信息）
        tree_obj.treedata = self._tree_to_treedata(phylo_tree)
        
        # 根据参数决定返回内容
        if only_return_tree:
            return tree_obj
        else:
            self.tree = tree_obj
            return self

    def _tree_to_treedata(self, tree: Tree) -> Dict:
        treedata = {
            'nodes': [],
            'edges': [],
            'root': None,
            'tips': []
        }
    
        # 1. 先获取所有节点的列表（固定顺序），并建立节点到索引的映射
        all_nodes = list(tree.find_clades())  # 转换为列表，固定迭代顺序
        node_to_idx = {node: idx for idx, node in enumerate(all_nodes)}  # 节点对象 -> 索引
    
        # 2. 收集所有节点信息
        for idx, node in enumerate(all_nodes):
            node_data = {
                'id': idx,
                'name': node.name,  # 叶节点有name（细胞ID），内部节点可能为None
                'branch_length': node.branch_length,
                'is_terminal': node.is_terminal()
            }
            treedata['nodes'].append(node_data)
        
            # 记录根节点
            if node is tree.root:
                treedata['root'] = idx
        
            # 记录末端节点（细胞）
            if node.is_terminal():
                treedata['tips'].append(idx)
    
        # 3. 收集分支信息（基于节点对象映射索引，避免依赖name）
        for parent in all_nodes:
            if not parent.is_terminal():  # 只处理内部节点（非叶节点）
                parent_idx = node_to_idx[parent]
                for child in parent.clades:  # 遍历子节点
                    child_idx = node_to_idx[child]
                    treedata['edges'].append({
                        'parent': parent_idx,
                        'child': child_idx,
                        'length': child.branch_length
                    })
    
        return treedata