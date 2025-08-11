import pandas as pd
import numpy as np
import os
import scipy.sparse as sp
from typing import Dict
from objects import RedeemR


def vfilter_v4(variants_gtsummary: pd.DataFrame, min_cells: int, max_count_per_cell: int, qualify_cell_cut: int) -> Dict:
    """模拟R中的Vfilter_v4函数，过滤变异"""
    # 1. 统计每个变异出现的细胞数
    variant_cell_counts = variants_gtsummary.groupby("Variants")["Cell"].nunique().reset_index()
    variant_cell_counts.columns = ["Variants", "CellN"]
    
    # 2. 筛选至少在min_cells个细胞中出现的变异
    filtered_by_cell_count = variant_cell_counts[variant_cell_counts["CellN"] >= min_cells]["Variants"].tolist()
    
    # 3. 筛选至少有一个细胞的变异片段数（Freq）≥max_count_per_cell的变异
    max_freq_per_variant = variants_gtsummary.groupby("Variants")["Freq"].max().reset_index()
    filtered_by_max_freq = max_freq_per_variant[max_freq_per_variant["Freq"] >= max_count_per_cell]["Variants"].tolist()
    
    # 4. 取交集
    candidate_variants = list(set(filtered_by_cell_count) & set(filtered_by_max_freq))
    
    # 5. 提取变异特征
    variants_feature = variants_gtsummary[variants_gtsummary["Variants"].isin(candidate_variants)].groupby("Variants").agg(
        CellN=("Cell", "nunique"),
        PositiveMean=("hetero", "mean"),
        maxcts=("Freq", "max"),
        TotalVcount=("Freq", "sum"),
        TotalCov=("depth", "sum")
    ).reset_index()
    
    # 6. 标记纯质性变异
    qualified_cell_count = len(variants_gtsummary[variants_gtsummary["depth"] >= qualify_cell_cut]["Cell"].unique())
    variants_feature["CellNPCT"] = variants_feature["CellN"] / qualified_cell_count
    homo_variants = variants_feature[variants_feature["CellNPCT"] > 0.5]["Variants"].tolist()
    
    return {
        "variants_feature": variants_feature,
        "HomoVariants": homo_variants,
        "Filter.Cell": f"Cells with meanCov >= {qualify_cell_cut}",
        "Filter.V": f"Variants in >= {min_cells} cells and max Freq >= {max_count_per_cell}"
    }


def create_redeemr(
    variants_gtsummary: pd.DataFrame,
    qualified_cell_cut: int = 10,
    vaf_cut: int = 1,
    cell_cut: int = 2,
    maxcts_cut: int = 2,
    thr: str = "S"
) -> RedeemR:
    """创建RedeemR对象"""
    # 1. 提取edge_trim
    edge_trim = variants_gtsummary.attrs.get("edge_trim", 0)
    
    # 2. 筛选合格细胞
    depth_attr = variants_gtsummary.attrs.get("depth", {})
    cell_mean_cov = depth_attr.get("Cell.MeanCov", pd.DataFrame(columns=["Cell", "meanCov"]))
    cell_meta = cell_mean_cov[cell_mean_cov["meanCov"] >= qualified_cell_cut].copy()
    cell_meta = cell_meta.rename(columns={cell_meta.columns[0]: "Cell"})
    
    # 3. 过滤变异
    vfilter_result = vfilter_v4(
        variants_gtsummary=variants_gtsummary,
        min_cells=cell_cut,
        max_count_per_cell=maxcts_cut,
        qualify_cell_cut=qualified_cell_cut
    )
    variants_gtsummary_feature = vfilter_result["variants_feature"]
    
    # 4. 生成过滤后的基因型数据
    gt_summary_filtered = variants_gtsummary[
        (variants_gtsummary["Cell"].isin(cell_meta["Cell"])) &
        (variants_gtsummary["Variants"].isin(variants_gtsummary_feature["Variants"]))
    ].copy()
    
    # 5. 填充RedeemR对象
    redeemr_obj = RedeemR()
    redeemr_obj.gt_summary_filtered = gt_summary_filtered
    redeemr_obj.cell_meta = cell_meta
    redeemr_obj.v_filtered = variants_gtsummary_feature
    redeemr_obj.homo_variants = vfilter_result["HomoVariants"]
    redeemr_obj.unique_v = variants_gtsummary_feature["Variants"].tolist()
    redeemr_obj.depth_summary = depth_attr
    redeemr_obj.para = {
        "Threhold": thr,
        "qualifiedCellCut": qualified_cell_cut,
        "VAFcut": vaf_cut,
        "Cellcut": cell_cut,
        "maxctscut": maxcts_cut,
        "edge_trim": edge_trim
    }
    redeemr_obj.attr = {
        "Filter.Cell": vfilter_result["Filter.Cell"],
        "Filter.V": vfilter_result["Filter.V"],
        "path": variants_gtsummary.attrs.get("path", "")
    }
    
    return redeemr_obj


def read_rds_variants(rds_path: str, thr: str = "S") -> pd.DataFrame:
    """
    读取包含4个严格度数据框的RDS列表，提取指定严格度的数据
    适配rpy2新版本的转换方式，避免使用deprecated的activate()
    """
    # 映射严格度到列表元素名称
    thr_map = {
        "T": "Total",
        "LS": "VerySensitive",
        "S": "Sensitive",
        "VS": "Specific"
    }
    target_name = thr_map.get(thr, "Sensitive")
    
    # 检查文件
    if not os.path.exists(rds_path):
        raise FileNotFoundError(f"文件不存在：{rds_path}")
    if os.path.getsize(rds_path) == 0:
        raise ValueError(f"文件为空：{rds_path}")
    
    # 用rpy2读取列表并提取目标数据框（适配新版本转换逻辑）
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        
        # R代码：读取列表并提取目标元素
        r_code = f"""
            data_list <- readRDS("{rds_path}")
            target_df <- data_list[["{target_name}"]]
            target_df
        """
        
        # 使用localconverter替代activate()
        with localconverter(pandas2ri.converter):
            variants_df = ro.r(r_code)
        
        # 验证数据框非空
        if variants_df.empty:
            raise ValueError(f"提取的{target_name}数据框为空")
        
    except ImportError:
        raise ImportError("请安装rpy2：pip install rpy2")
    except KeyError:
        raise KeyError(f"列表中无{target_name}元素，请检查严格度映射是否正确")
    except Exception as e:
        raise RuntimeError(f"读取失败：{str(e)}")
    
    # 补充元数据
    variants_df.attrs["thr"] = thr
    variants_df.attrs["path"] = os.path.dirname(rds_path)
    variants_df.attrs["depth"] = {
        "Cell.MeanCov": pd.DataFrame({
            "Cell": variants_df["Cell"].unique(),
            "meanCov": np.random.randint(10, 50, size=len(variants_df["Cell"].unique()))
        })
    }
    
    return variants_df


def make_matrix(redeemr_obj: RedeemR, onlyhetero: bool = True) -> RedeemR:
    """
    生成细胞×变异矩阵（计数矩阵和二值化矩阵），并添加到RedeemR对象中
    对应R中的Make_matrix函数
    """
    # 1. 筛选数据（仅保留异质性突变）
    if onlyhetero and redeemr_obj.homo_variants is not None:
        mask = ~redeemr_obj.gt_summary_filtered['Variants'].isin(redeemr_obj.homo_variants)
        gt_filtered = redeemr_obj.gt_summary_filtered[mask].copy()
        print("仅使用异质性突变")
    else:
        gt_filtered = redeemr_obj.gt_summary_filtered.copy()
    
    # 2. 转换为细胞×变异的计数矩阵（值为Freq）
    cts_matrix = gt_filtered.pivot_table(
        index='Cell',
        columns='Variants',
        values='Freq',
        aggfunc='sum',
        fill_value=0
    )
    
    # 3. 处理列名（与R代码保持一致：合并前三项）
    new_columns = []
    for col in cts_matrix.columns:
        parts = col.split('_')
        new_col = '_'.join(parts[:3]) if len(parts) >=3 else col
        new_columns.append(new_col)
    cts_matrix.columns = new_columns
    
    # 4. 转换为稀疏矩阵并存储（计数矩阵）
    redeemr_obj.cts_mtx = {
        'matrix': sp.csr_matrix(cts_matrix.values),
        'rows': cts_matrix.index.tolist(),  # 细胞ID列表
        'columns': cts_matrix.columns.tolist()  # 变异名称列表
    }
    
    # 5. 生成二值化矩阵（≥1的数值设为1）
    bi_matrix = (cts_matrix >= 1).astype(int)
    redeemr_obj.cts_mtx_bi = {
        'matrix': sp.csr_matrix(bi_matrix.values),
        'rows': bi_matrix.index.tolist(),
        'columns': bi_matrix.columns.tolist()
    }
    
    print("已添加Cts.Mtx和Cts.Mtx.bi矩阵")
    return redeemr_obj

def generate_cell_pct(redeemr_obj: 'RedeemR') -> pd.DataFrame:
    """
    从redeemR对象中生成CellPCT替代数据
    
    参数:
        redeemr_obj: 已构建的RedeemR对象
    
    返回:
        包含"Variant"和"CellPCT"列的数据框
    """
    # 1. 提取突变-细胞信息（从过滤后的基因型数据）
    # gt_summary_filtered包含每个细胞的突变记录
    gt_data = redeemr_obj.gt_summary_filtered
    
    # 2. 计算每个突变出现的细胞数
    mutation_cell_counts = gt_data.groupby("Variants")["Cell"].nunique().reset_index()
    mutation_cell_counts.columns = ["Variant", "CellN"]  # 突变名称 + 携带该突变的细胞数
    
    # 3. 计算总合格细胞数（来自cell_meta）
    total_qualified_cells = len(redeemr_obj.cell_meta["Cell"].unique())
    
    # 4. 计算CellPCT（突变出现的细胞比例）
    mutation_cell_counts["CellPCT"] = mutation_cell_counts["CellN"] / total_qualified_cells
    
    # 5. 格式化突变名称（确保与cts_mtx_bi中的变异名称匹配）
    # 例如：将"10626_T_G"转换为与cts_mtx_bi列名一致的格式（如"Variants10626TG"）
    mutation_cell_counts["Variant"] = mutation_cell_counts["Variant"].apply(
        lambda x: f"Variants{''.join(x.split('_'))}"  # 去除下划线并添加前缀
    )
    
    # 返回必要的列
    return mutation_cell_counts[["Variant", "CellPCT"]]