import os
import scanpy as sc
import snapatac2 as snap
import muon as mu
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy import stats

# 设置中文字体
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]

# 创建日志目录和输出目录
output_dir = "/data/project/r2py/young2-result"
log_dir = os.path.join(output_dir, "logs")
figures_dir = os.path.join(output_dir, "figures")
os.makedirs(log_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

# 日志记录函数
def log(message):
    timestamp = pd.Timestamp.now().strftime('%H:%M:%S')
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    
    # 写入日志文件
    log_file = os.path.join(log_dir, "analysis.log")
    with open(log_file, "a") as f:
        f.write(log_msg + "\n")

log("开始RNA-seq与ATAC-seq联合分析...")

# =================================
# 第一部分：加载并预处理单组学数据
# =================================
log("加载单组学数据...")

try:
    # 加载RNA数据（假设已完成基础预处理）
    rna_adata = sc.read_h5ad("/data/project/r2py/dataset/Young2.HSC_processed.h5ad")
    log(f"RNA数据加载完成: {rna_adata.shape}")
    
    # 加载ATAC数据（假设已完成基础预处理）
    atac_adata = snap.read("/data/project/r2py/result/initial_data.h5ad")
    log(f"ATAC数据加载完成: {atac_adata.shape}")
    
    # 检查关键列是否存在
    required_cols = ['chr', 'start', 'end']
    missing_cols = []
    
    for col in required_cols:
        try:
            _ = atac_adata.var[col]
        except KeyError:
            missing_cols.append(col)
    
    if missing_cols:
        log(f"警告：ATAC数据缺少染色体位置信息：{', '.join(missing_cols)}")
        raise ValueError("缺少必要的染色体位置信息")
    
except Exception as e:
    log(f"数据加载失败: {e}")
    exit(1)

# =================================
# 第二部分：数据对齐与多组学对象创建
# =================================
log("开始数据对齐...")

# 细胞ID对齐
common_cells = rna_adata.obs_names.intersection(atac_adata.obs_names)
common_cells = list(common_cells)

if len(common_cells) == 0:
    log("错误：RNA和ATAC数据中没有共同的细胞ID")
    exit(1)

# RNA数据切片
rna_adata = rna_adata[common_cells, :]

# ATAC数据切片
atac_adata = atac_adata.subset(
    obs_indices=common_cells,
    inplace=False
)
log(f"匹配后共保留 {len(common_cells)} 个细胞")

# 创建多组学对象
log("创建多组学对象...")
mdata = mu.MuData({"rna": rna_adata, "atac": atac_adata})

# =================================
# 第三部分：ATAC数据基因注释与活性计算
# =================================
log("开始ATAC数据基因注释...")

try:
    # 计算基因活性矩阵
    gene_activity_adata = snap.pp.make_gene_matrix(
        mdata['atac'],
        gene_anno="/data/project/r2py/dataset/gencode_v41_GRCh38.gff3.gz",
        use_x=True,
        inplace=False,
    )
    
    # 将基因活性矩阵作为新模态添加到多组学对象
    mdata.mod["gene_activity"] = gene_activity_adata
    log("基因活性矩阵计算完成，已作为新模态添加到多组学对象")
    
except Exception as e:
    log(f"基因注释或活性计算失败: {e}")
    exit(1)

# =================================
# 第四部分：联合降维与聚类分析
# =================================
log("开始联合降维与聚类分析...")

try:
    # 提取RNA表达矩阵和ATAC基因活性矩阵
    rna_expr = mdata['rna'].X
    gene_activity = mdata['gene_activity'].X
    
    # 转换为密集矩阵（如果是稀疏矩阵）
    if hasattr(rna_expr, 'toarray'):
        rna_expr = rna_expr.toarray()
    if hasattr(gene_activity, 'toarray'):
        gene_activity = gene_activity.toarray()
    
    # 标准化处理
    rna_norm = sc.pp.scale(rna_expr, copy=True)
    activity_norm = sc.pp.scale(gene_activity, copy=True)
    
    # 特征匹配：找到两个数据集中共同的基因
    common_genes = list(set(mdata['rna'].var_names).intersection(set(mdata['gene_activity'].var_names)))
    log(f"找到 {len(common_genes)} 个共同基因用于联合分析")
    
    if len(common_genes) == 0:
        log("错误：RNA表达矩阵和基因活性矩阵中没有共同基因")
        exit(1)
    
    # 提取共同基因的表达和活性
    rna_common_idx = [mdata['rna'].var_names.get_loc(gene) for gene in common_genes]
    activity_common_idx = [mdata['gene_activity'].var_names.get_loc(gene) for gene in common_genes]
    
    rna_common = rna_norm[:, rna_common_idx]
    activity_common = activity_norm[:, activity_common_idx]
    
    # 合并特征并进行PCA降维
    combined_features = np.hstack([rna_common, activity_common])
    pca = PCA(n_components=50)
    mdata['rna'].obsm['X_joint_pca'] = pca.fit_transform(combined_features)
    
    # 计算PCA方差解释率
    variance_ratio = pca.explained_variance_ratio_
    log(f"PCA降维完成，前50个主成分方差解释率: {sum(variance_ratio):.4f}")
    
    # 可视化PCA方差解释率
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, len(variance_ratio) + 1), variance_ratio, 'o-')
    plt.xlabel('主成分')
    plt.ylabel('方差解释率')
    plt.title('PCA方差解释率')
    plt.savefig(os.path.join(figures_dir, 'pca_variance_ratio.png'))
    plt.close()
    
    # 基于联合PCA结果进行UMAP和聚类
    sc.pp.neighbors(mdata['rna'], use_rep='X_joint_pca', n_neighbors=15)
    sc.tl.umap(mdata['rna'], min_dist=0.3)
    
    # 尝试不同分辨率的聚类
    resolutions = [0.3, 0.5, 0.8]
    for res in resolutions:
        key = f'joint_leiden_r{res}'
        sc.tl.leiden(mdata['rna'], resolution=res, key_added=key)
    
    # 可视化联合聚类结果
    fig = sc.pl.umap(
        mdata['rna'],
        color=[f'joint_leiden_r{r}' for r in resolutions],
        ncols=3,
        title=[f'分辨率: {r}' for r in resolutions],
        return_fig=True
    )
    save_path = os.path.join(figures_dir, 'joint_clustering.png')
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    log(f"联合聚类完成，尝试了{len(resolutions)}种不同分辨率，聚类图已保存至 {save_path}")
    
except Exception as e:
    log(f"联合降维或聚类失败: {e}")
    exit(1)

# =================================
# 第五部分：关联分析与可视化
# =================================
log("开始关联分析...")

def correlate_rna_atac(mdata, gene_list):
    results = []
    for gene in gene_list:
        try:
            # 获取RNA表达量
            if gene not in mdata['rna'].var_names:
                log(f"警告: RNA表达矩阵中未找到基因 {gene}")
                continue
            rna_expr = mdata['rna'].obs_vector(gene)
            
            # 获取基因活性
            if gene not in mdata['gene_activity'].var_names:
                log(f"警告: 基因活性矩阵中未找到基因 {gene}")
                continue
            atac_signal = mdata['gene_activity'].obs_vector(gene)
            
            # 计算相关性
            corr, pval = stats.spearmanr(rna_expr, atac_signal)
            results.append({"gene": gene, "correlation": corr, "pvalue": pval})
        except Exception as e:
            log(f"计算基因 {gene} 相关性时出错: {e}")
    
    return pd.DataFrame(results)

try:
    # 选择高变基因进行分析
    top_variable_genes = mdata['rna'].var.sort_values('dispersions_norm', ascending=False).head(10).index
    corr_df = correlate_rna_atac(mdata, top_variable_genes)
    
    # 保存结果
    corr_df.to_csv(os.path.join(output_dir, 'gene_atac_correlation.csv'), index=False)
    log(f"相关性分析完成，已保存{len(corr_df)}个基因的结果")
    
    # 可视化高相关性基因
    if not corr_df.empty:
        top_genes = corr_df.sort_values('correlation', ascending=False).head(3)['gene'].tolist()
        for gene in top_genes:
            try:
                # 确保UMAP已计算
                if 'X_umap' not in mdata['rna'].obsm:
                    sc.tl.umap(mdata['rna'])
                
                # 准备基因活性数据用于可视化
                mdata['rna'].obs[f"gene_activity_{gene}"] = mdata['gene_activity'].obs_vector(gene)
                
                # 保存UMAP图像
                fig_path = os.path.join(figures_dir, f"{gene}_expression_activity")
                fig_abs_path = os.path.abspath(fig_path)
                sc.pl.umap(
                    mdata['rna'],
                    color=[gene, f"gene_activity_{gene}"],
                    title=[f"{gene}表达量", f"{gene}基因活性"],
                    ncols=2,
                    save=fig_abs_path  # scanpy会自动添加.png后缀
                )
                log(f"成功保存 {gene} 的UMAP可视化结果")
            except Exception as e:
                log(f"可视化 {gene} 时出错: {e}")
    else:
        log("警告: 相关性分析结果为空，未进行可视化")
    
except Exception as e:
    log(f"关联分析失败: {e}")

# =================================
# 第六部分：保存结果
# =================================
log("保存分析结果...")

try:
    # 为每个模态的var添加唯一前缀
    mdata['rna'].var_names = [f"rna_{gene}" for gene in mdata['rna'].var_names]
    mdata['atac'].var_names = [f"atac_{gene}" for gene in mdata['atac'].var_names]
    mdata['gene_activity'].var_names = [f"activity_{gene}" for gene in mdata['gene_activity'].var_names]
    
    # 保存多组学分析结果
    mdata.write(os.path.join(output_dir, "rna_atac_combined.h5mu"))
    log(f"多组学分析结果已保存至: {os.path.join(output_dir, 'rna_atac_combined.h5mu')}")
    
    # 保存单组学结果（带联合分析标签）
    mdata['rna'].write(os.path.join(output_dir, "rna_combined.h5ad"))
    mdata['atac'].write(os.path.join(output_dir, "atac_combined.h5ad"))
    log("单组学结果已保存")
    
    log("RNA-seq与ATAC-seq联合分析全部完成!")
    
except Exception as e:
    log(f"结果保存失败: {e}")