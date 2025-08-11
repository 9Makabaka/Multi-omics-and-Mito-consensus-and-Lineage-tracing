import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]  # 支持中文的字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

result_dir = "/data/project/r2py/young2RNA-result"
os.makedirs(result_dir, exist_ok=True)  # 若目录不存在则创建，存在则不报错
print(f"所有结果将保存至: {result_dir}")

"""
Step 1: 数据加载和检查
"""
# 加载H5AD文件
adata = sc.read_h5ad("/data/dataset/Young2.HSC.h5ad")

"""
Step 2: 质量控制计算
"""
# 标记特殊基因集
print("\n=== 标记质量控制基因 ===")
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))  # 线粒体基因
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))  # 核糖体基因
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")  # 血红蛋白基因

# 确保QC指标计算
print("\n=== 计算QC指标 ===")
if 'n_genes' not in adata.obs:
    # 计算每个细胞中表达的基因数（表达>0）
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # 计算总表达量
    adata.obs["total_expression"] = np.array(adata.X.sum(axis=1)).flatten()
    
    # 计算线粒体基因百分比
    if adata.var["mt"].sum() > 0:
        mt_expression = np.array(adata.X[:, adata.var["mt"]].sum(axis=1)).flatten()
        adata.obs["pct_mt"] = mt_expression / adata.obs["total_expression"] * 100
    else:
        adata.obs["pct_mt"] = 0
        print("警告: 未检测到线粒体基因!")

# 打印确认
print("已计算的QC指标:", list(adata.obs.columns[-3:]))

"""
Step 3: 质量控制可视化
"""
def plot_qc_metrics(adata, save_path):
    # 检查必要指标是否存在
    required_metrics = ["n_genes", "total_expression", "pct_mt"]
    missing = [m for m in required_metrics if m not in adata.obs]
    if missing:
        raise ValueError(f"缺少必要的QC指标: {missing}")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    metrics = required_metrics
    titles = ["表达基因数", "总表达量", "线粒体基因百分比"]
    
    for ax, metric, title in zip(axes, metrics, titles):
        data = adata.obs[metric]
        
        # 绘制简化版小提琴图
        parts = ax.violinplot(
            [data],
            showmeans=False,
            showmedians=True,
            widths=0.7
        )
        
        # 设置颜色
        for pc in parts['bodies']:
            pc.set_facecolor('#1f77b4')
            pc.set_edgecolor('black')
            pc.set_alpha(0.7)
        
        parts['cmedians'].set_color('red')
        parts['cmins'].set_color('black')
        parts['cmaxes'].set_color('black')
        parts['cbars'].set_color('black')
        
        # 添加标题和标签
        ax.set_title(title)
        ax.set_ylabel(metric)
        ax.set_xticks([])
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')  # 保存图片
    plt.close()  # 关闭画布释放内存
    print(f"QC指标图已保存至: {save_path}")


print("\n=== 可视化QC指标 ===")
qc_fig_path = os.path.join(result_dir, "qc_metrics_violin.png")
plot_qc_metrics(adata, qc_fig_path)

"""
Step 3: 细胞和基因过滤 (调整方法)
"""
print("\n=== 过滤低质量细胞和基因 ===")
# 保存当前数据到raw
adata.raw = adata

# 细胞过滤 (基于表达基因数)
min_genes = 100
max_genes = np.percentile(adata.obs["n_genes"], 99.5)  # 过滤极端高表达细胞
sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["n_genes"] < max_genes, :]

# 基因过滤 (在至少3个细胞中表达)
sc.pp.filter_genes(adata, min_cells=3)

print(f"过滤后数据集形状: {adata.shape}")

"""
Step 4: 高变基因选择 (使用标准化数据的方法)
"""
print("\n=== 高变基因选择 ===")
if 'highly_variable' not in adata.var:
    # 方法1: 尝试Seurat的vst方法
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=2000,
            flavor="seurat_v3",
            layer=None  # 使用标准化数据
        )
    except Exception as e:
        print(f"Seurat方法失败: {e}")
        # 方法2: 使用基于分散度的方法
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=2000,
            flavor="cell_ranger"
        )
        
# 可视化高变基因并保存
hv_fig_path = os.path.join(result_dir, "highly_variable_genes.png")
sc.pl.highly_variable_genes(
    adata,
    save=False,
    show=False  # 不自动显示
)
plt.savefig(hv_fig_path, dpi=300, bbox_inches='tight')  # 手动保存
plt.close()  # 关闭画布
# 检查文件是否存在
if os.path.exists(hv_fig_path):
    print(f"高变基因可视化图已保存至: {hv_fig_path}")
else:
    print(f"警告：高变基因可视化图保存失败！")

"""
Step 5: 数据缩放和PCA
"""
print("\n=== 数据缩放和PCA ===")
# 如果数据尚未缩放
if 'scale' not in adata.layers:
    sc.pp.scale(adata, max_value=10)
    adata.layers["scale"] = adata.X.copy()

# PCA分析
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata, use_highly_variable=True, svd_solver='arpack')
    # 保存PCA方差比例图
    pca_var_path = os.path.join(result_dir, "pca_variance_ratio.png")
    sc.pl.pca_variance_ratio(
        adata, 
        n_pcs=50, 
        log=True,
        save=False,
        show=False
    )
    plt.savefig(pca_var_path, dpi=300, bbox_inches='tight')  # 手动保存
    plt.close()  # 关闭画布
    # 检查文件是否存在
    if os.path.exists(pca_var_path):
        print(f"PCA方差比例图已保存至: {pca_var_path}")
    else:
        print(f"警告：PCA方差比例图保存失败！")

"""
Step 6: 邻域图和UMAP可视化
"""
print("\n=== 邻域图和UMAP ===")
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)

if 'X_umap' not in adata.obsm:
    sc.tl.umap(adata)
    # 保存UMAP图（按QC指标着色）
    umap_qc_path = os.path.join(result_dir, "umap_qc_metrics.png")
    sc.pl.umap(
        adata, 
        color=["n_genes", "total_expression", "pct_mt"], 
        ncols=3,
        save=False,
        show=False
    )
    # 调整路径
    plt.savefig(umap_qc_path, dpi=300, bbox_inches='tight')  # 手动保存
    plt.close()  # 关闭画布
    # 检查文件是否存在
    if os.path.exists(umap_qc_path):
        print(f"UMAP质量指标图已保存至: {umap_qc_path}")
    else:
        print(f"警告：UMAP质量指标图保存失败！")

"""
Step 7: 聚类分析
"""
print("\n=== 聚类分析 ===")
if 'leiden' not in adata.obs:
    # 尝试不同分辨率
    for res in [0.5, 0.8, 1.0]:
        key = f"leiden_r{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        # 保存当前分辨率的聚类UMAP图
        cluster_fig_path = os.path.join(result_dir, f"umap_leiden_r{res}.png")
        sc.pl.umap(
            adata, 
            color=key, 
            title=f"Leiden r={res}", 
            legend_loc='on data',
            show=False,  # 不显示图片
            save=False   # 禁用默认save参数
        )
        plt.savefig(cluster_fig_path, dpi=300, bbox_inches='tight')  # 手动保存
        plt.close()  # 关闭画布
        # 检查文件是否存在
        if os.path.exists(cluster_fig_path):
            print(f"Leiden聚类图(r={res})已保存至: {cluster_fig_path}")
        else:
            print(f"警告：Leiden聚类图(r={res})保存失败！")

"""
Step 8: 差异表达分析
"""
print("\n=== 差异表达分析 ===")
if 'rank_genes_groups' not in adata.uns:
    # 使用中等分辨率聚类结果
    cluster_key = "leiden_r0.5" if "leiden_r0.5" in adata.obs else "leiden"
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')
    
    # 保存marker基因dotplot
    dotplot_path = os.path.join(result_dir, "marker_genes_dotplot.png")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_key,
        n_genes=5,
        standard_scale='var',
        swap_axes=True,
        save=False,
        show=False
    )
    # 调整路径
    plt.savefig(dotplot_path, dpi=300, bbox_inches='tight')  # 手动保存
    plt.close()  # 关闭画布
    # 检查文件是否存在
    if os.path.exists(dotplot_path):
        print(f"Marker基因dotplot已保存至: {dotplot_path}")
    else:
        print(f"警告：Marker基因dotplot保存失败！")

"""
Step 9: 保存结果
"""
output_file = os.path.join(result_dir, "Young2.HSC_processed.h5ad")
adata.write_h5ad(output_file)
print(f"\n分析完成! 结果已保存至: {output_file}")