# R2py

# 一、`young2_RNA.py`

[young2_RNA.py](young2_RNA.py)

## **1.输入文件**

- 输入数据：`/data/dataset/Young2.HSC.h5ad`这是一个 H5AD 格式文件，存储单细胞 RNA 测序（scRNA-seq）原始数据，包含基因表达矩阵、细胞元数据（obs）和基因元数据（var）等信息。

## **2.输出文件**

1. **`qc_metrics_violin.png`**
    - 内容：包含 3 个小提琴图，分别展示所有细胞的 “表达基因数”“总表达量”“线粒体基因百分比” 的分布。
    - 作用：直观呈现细胞质量的基础分布，帮助判断数据整体质量（如是否存在大量低质量细胞）。
2. **`highly_variable_genes.png`**
    - 内容：高变基因筛选结果的可视化图，通常以 “平均表达量” 为 x 轴、“分散度（或标准化分散度）” 为 y 轴，标记出被筛选为 “高变基因” 的点（通常用红色标注）。
    - 作用：展示高变基因的筛选依据，验证筛选结果的合理性（高变基因是后续降维、聚类的核心特征）。
3. **`pca_variance_ratio.png`**
    - 内容：前 50 个主成分（PC）的方差解释比例折线图（y 轴为方差比例，x 轴为主成分编号），通常采用对数坐标（log=True）。
    - 作用：判断 PCA 降维的有效性，通过 “肘部法则” 确定后续分析中应保留的主成分数量（通常取方差解释比例趋于平缓的拐点前的 PCs）。
4. **`umap_qc_metrics.png`**
    - 内容：UMAP 降维后的 2D 散点图（3 个子图并排），分别按 “表达基因数（n_genes）”“总表达量（total_expression）”“线粒体基因百分比（pct_mt）” 着色。
    - 作用：展示细胞在低维空间的分布与质量指标的关联（如低质量细胞是否聚集在特定区域）。
5. **`umap_leiden_r0.5.png`、`umap_leiden_r0.8.png`、`umap_leiden_r1.0.png`**
    - 内容：UMAP 降维后的散点图，按对应分辨率（0.5、0.8、1.0）的 Leiden 聚类结果着色，每个颜色代表一个细胞亚群，图例标注亚群编号。
    - 作用：展示不同分辨率下细胞亚群的划分结果（分辨率越高，聚类越精细，亚群数量越多），用于选择合适的亚群划分方案。
6. **`marker_genes_dotplot.png`**
    - 内容：点图（dot plot），横轴为细胞亚群（基于分辨率 0.5 的 Leiden 聚类结果），纵轴为每个亚群的前 5 个 marker 基因（差异表达基因）。点的大小表示该基因在亚群中表达的细胞比例，点的颜色表示该基因在亚群中的平均表达量（通常经标准化）。
    - 作用：直观展示各细胞亚群的特征基因，用于解释亚群的生物学意义（如判断亚群对应的细胞类型）。
7. **`Young2.HSC_processed.h5ad`**
    - 内容：整合了所有分析结果的数据容器，包括：
        - 过滤后的基因表达矩阵（`adata.X`）；
        - 细胞元数据（`adata.obs`）：包含 QC 指标（n_genes、total_expression、pct_mt）、不同分辨率的 Leiden 聚类标签等；
        - 基因元数据（`adata.var`）：包含基因是否为高变基因（`highly_variable`）等信息；
        - 降维结果（`adata.obsm`）：PCA 坐标（`X_pca`）、UMAP 坐标（`X_umap`）；
        - 差异表达分析结果（`adata.uns['rank_genes_groups']`）：每个亚群的差异基因排名、log2 倍数变化、统计显著性等。
    - 作用：保存完整的分析中间结果和最终结果，便于后续深入分析（如细胞类型注释、通路富集分析等）或结果复现。

## **3.处理步骤**

基于`scanpy`库实现了 scRNA-seq 数据的完整预处理和基础分析流程。

### **3.1数据加载和初始化**

- 加载 H5AD 格式的原始单细胞数据到`adata`对象（scanpy 的核心数据结构，整合表达矩阵、元数据等）。

### **3.2质量控制指标计算**

- **标记特征基因**：识别线粒体基因（`mt`，名称以`MT-`或`mt-`开头）、核糖体基因（`ribo`，名称以`RPS`或`RPL`开头）、血红蛋白基因（`hb`，名称含`HB`且不含`P`），用于后续 QC。
- **计算 QC 指标**：
    - `n_genes`：每个细胞中表达的基因数（表达量 > 0 的基因）；
    - `total_expression`：每个细胞的总表达量（所有基因表达量之和）；
    - `pct_mt`：线粒体基因表达量占细胞总表达量的百分比（用于判断细胞质量，高比例可能提示细胞受损）。

```python
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
```

### **3. QC 指标可视化**

- 通过`plot_qc_metrics`函数绘制小提琴图，展示`n_genes`、`total_expression`、`pct_mt`的分布，直观判断数据质量（如是否存在低质量细胞）。

### **4. 细胞和基因过滤**

- **细胞过滤**：保留表达基因数在 100 到 99.5 百分位之间的细胞（过滤基因数过少的低质量细胞和基因数异常高的可能 doublet 细胞）。
- **基因过滤**：保留至少在 3 个细胞中表达的基因（过滤低丰度、可能不可靠的基因）。
- 保存原始数据到`adata.raw`，用于后续差异分析时回溯原始表达量。

```python
# 细胞过滤 (基于表达基因数)
min_genes = 100
max_genes = np.percentile(adata.obs["n_genes"], 99.5)  # 过滤极端高表达细胞
sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["n_genes"] < max_genes, :]

# 基因过滤 (在至少3个细胞中表达)
sc.pp.filter_genes(adata, min_cells=3)

### 过滤后数据集形状: (3683, 3000)
```

### **5. 高变基因选择**

- 识别 2000 个高变基因（在不同细胞中表达差异显著的基因，是后续分析的核心特征）：
    - 优先使用 Seurat 的`vst`方法；若失败， fallback 到基于分散度的`cell_ranger`方法。
- 可视化高变基因的分散度（表达量与分散度的关系），验证高变基因选择合理性。

```python
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
```

### **6. 数据缩放和 PCA 降维**

- **数据缩放**：对基因表达量进行标准化（均值为 0，方差为 1），避免高表达基因主导分析，结果保存到`scale`层。
- **PCA 分析**：基于高变基因进行主成分分析（PCA），提取前 30 个主成分（PCs），用于后续降维和聚类；绘制 PCA 方差比例图，判断主成分的解释力。

```python
# 如果数据尚未缩放
if 'scale' not in adata.layers:
    sc.pp.scale(adata, max_value=10)
    adata.layers["scale"] = adata.X.copy()

# PCA分析
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata, use_highly_variable=True, svd_solver='arpack')
```

### **7. 邻域图构建和 UMAP 可视化**

- **邻域图**：基于 PCA 结果构建细胞间的近邻关系图（使用 30 个 PCs 和 20 个近邻），反映细胞间的相似性。
- **UMAP 降维**：通过 UMAP 算法将高维 PCA 结果降维到 2D 空间，实现细胞群体的可视化；按`n_genes`、`total_expression`、`pct_mt`着色，观察质量指标与细胞分布的关系。

```python
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)

if 'X_umap' not in adata.obsm:
    sc.tl.umap(adata)
```

### **8. 聚类分析**

- 使用 Leiden 算法对细胞进行聚类（基于邻域图），尝试 3 种分辨率（0.5、0.8、1.0），获得不同精细度的细胞亚群；在 UMAP 图上可视化聚类结果，展示不同亚群的分布。

```python
for res in [0.5, 0.8, 1.0]:
    key = f"leiden_r{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key)
```

### **9. 差异表达分析**

- 以中等分辨率（0.5）的聚类结果为分组，用 Wilcoxon 秩和检验计算每个亚群的差异表达基因（marker 基因）。
- 绘制_dotplot 展示每个亚群的前 5 个 marker 基因，直观反映基因在不同亚群中的表达特异性。

```python
cluster_key = "leiden_r0.5" if "leiden_r0.5" in adata.obs else "leiden"
sc.tl.rank_genes_groups(adata, groupby=cluster_key, method='wilcoxon')
```

### **10. 结果保存**

- 将处理后的`adata`对象（包含所有分析结果）写入 H5AD 文件，供后续深入分析或可视化使用。

## 4.输出结果

`qc_metrics_violin.png`

![qc_metrics_violin.png](qc_metrics_violin.png)

`highly_variable_genes.png`

![highly_variable_genes.png](highly_variable_genes.png)

`pca_variance_ratio.png`

![pca_variance_ratio.png](pca_variance_ratio.png)

`marker_genes_dotplot.png`

![marker_genes_dotplot.png](marker_genes_dotplot.png)

`umap_qc_metrics.png`

![umap_qc_metrics.png](umap_qc_metrics.png)

`umap_leiden_r0.5.png`

![umap_leiden_r0.5.png](umap_leiden_r0.5.png)

`umap_leiden_r0.8.png`

![umap_leiden_r0.8.png](umap_leiden_r0.8.png)

`umap_leiden_r1.0.png`

![umap_leiden_r1.0.png](umap_leiden_r1.0.png)

# 二、`young2_ATAC.py`

[young2_ATAC.py](young2_ATAC.py)

## **1.输入文件**

1. **`atac_counts.mtx`**
    - 格式：稀疏矩阵文件（Matrix Market 格式）
    - 内容：单细胞 ATAC-seq 的 peak 计数矩阵， rows 对应细胞，columns 对应 peak（染色质可及性区域），值为每个细胞中每个 peak 的测序计数。
2. **`peak_annotations.csv`**
    - 格式：CSV 文件
    - 内容：peak 的注释信息（如染色体位置、基因关联等），索引为 peak ID，列包含 peak 的属性（如`chr`、`start`、`end`等，用于定位 peak 在基因组上的位置）。
3. **`cell_metadata.csv`**
    - 格式：CSV 文件
    - 内容：细胞的元数据（如测序质量指标、样本来源等），索引为细胞 ID，列包含细胞的属性（如`n_counts`、`tsse`等，用于描述细胞的基础特征）。

## **2.输出文件**

1. **`initial_data.h5ad`**
    - 格式：AnnData 对象（单细胞分析的标准数据结构）
    - 内容：整合后的 ATAC-seq 数据，包含：
        - 计数矩阵（`adata.X`）、细胞元数据（`adata.obs`）、peak 注释（`adata.var`）；
        - 降维和聚类结果（如 UMAP 坐标、Leiden 聚类标签）。
2. **可视化图片**
    - `umap_initial.png`：UMAP 降维图，按细胞的基础指标（`n_counts`、`tsse`）着色，展示细胞分布与质量指标的关系；
    - `umap_leiden.png`：UMAP 降维图，按 Leiden 聚类结果着色，展示不同细胞簇的分布。

## **3.处理步骤**

**单细胞 ATAC-seq 数据的基础预处理与分析流程。**

### **3.1数据读取与预处理**

- **读取原始数据**：
    - 读取稀疏矩阵`atac_counts.mtx`并转置（确保 rows 为细胞、columns 为 peak）；
    - 读取 peak 注释（`peak_annotations.csv`）和细胞元数据（`cell_metadata.csv`）。
- **数据校验**：
    - 检查矩阵形状与元数据 / 注释的数量是否匹配（避免细胞 /peak 数量不一致）；
    - 处理元数据中的缺失值（数值列填充 0，字符串列填充`'nan'`）；
    - 转换元数据列类型（确保符合 AnnData 要求：数值型、布尔型或分类型，避免不支持的类型）。

### **3.2创建 AnnData 对象**

- 基于读取的数据构建`AnnData`对象（单细胞分析的标准容器），整合：
    - 计数矩阵（`adata.X`）：存储细胞 - peak 的计数信息；
    - 细胞元数据（`adata.obs`）：存储细胞的属性（如质量指标）；
    - peak 注释（`adata.var`）：存储 peak 的基因组位置等信息（若包含`chr`、`start`、`end`列，则单独标记为基因组坐标）。
    
    ```python
    AnnData object with n_obs x n_vars = 3703 x 214525 backed at '/data/project/r2py/young2ATAC-result/initial_data.h5ad'
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC', 'nCount_SCT', 'nFeature_SCT', 'SCT.weight', 'ATAC.weight', 'seurat_clusters', 'STD.CellType', 'STD_Cat', 'STD_Cat2', 'Sample', 'Filter', 'HLF', 'CRHBP', 'CD34', 'Rig.HSC', 'MitoCoverage', 'ClonalGroup'
        var: 'chr', 'start', 'end'
    AnnData对象创建完成，包含 3703 个细胞和 214525 个特征
    ```
    

### **3.3降维分析**

- **特征选择**：通过`snap.pp.select_features`选择 250,000 个高变 peak（聚焦于可及性差异显著的区域，减少噪声）；
- **谱嵌入（Spectral Embedding）**：一种降维方法，保留数据的全局结构，为后续可视化和聚类提供基础；
- **UMAP 降维**：将高维数据映射到 2D 空间，便于可视化细胞群体分布；
- **可视化**：生成`umap_initial.png`，按细胞的基础指标（如测序计数`n_counts`、TSSE 值）着色，观察细胞质量与分布的关系。

### **3.4聚类分析**

- **构建 K 近邻图**：基于降维后的细胞相似性构建图结构（细胞为节点，相似性为边）；
- **Leiden 聚类**：基于图结构的聚类算法，将相似的细胞聚为一簇（识别细胞亚群）；
- **可视化**：生成`umap_leiden.png`，按聚类标签着色，展示不同细胞簇的空间分布，并统计簇的数量。

```python
聚类完成，共识别出 7 个细胞簇
```

## 4.输出结果

`umap_initial.png`

![umap_initial.png](umap_initial.png)

`umap_leiden.png`

![umap_leiden.png](umap_leiden.png)

# 三、`young2_RNA-ATAC.py`

[young2_RNA-ATAC.py](young2_RNA-ATAC.py)

## **1.输入文件**

1. **RNA-seq 数据**
    - 文件路径：`/data/project/r2py/dataset/Young2.HSC_processed.h5ad`
    - 格式：AnnData 对象（.h5ad）
    - 说明：已完成基础预处理的 RNA-seq 数据，包含基因表达矩阵和细胞元数据
2. **ATAC-seq 数据**
    - 文件路径：`/data/project/r2py/result/initial_data.h5ad`
    - 格式：AnnData 对象（.h5ad）
    - 说明：已完成基础预处理的 ATAC-seq 数据，包含染色质可及性矩阵和细胞元数据
3. **基因组注释文件**
    - 文件路径：`/data/project/r2py/dataset/gencode_v41_GRCh38.gff3.gz`
    - 格式：GFF3 压缩文件
    - 说明：用于 ATAC 数据的基因注释，计算基因活性矩阵

## **2.输出文件**

1. **日志文件**
    - 文件路径：`/data/project/r2py/young2-result/logs/analysis.log`
    - 格式：文本文件
    - 说明：记录整个分析过程的详细日志，包括每个步骤的开始时间、完成情况和可能的错误信息
2. **可视化结果**
    - 保存目录：`/data/project/r2py/young2-result/figures/`
    - 包含以下图像文件：
        - `pca_variance_ratio.png`：PCA 降维的方差解释率图
        - `joint_clustering.png`：不同分辨率下的联合聚类 UMAP 图
        - `{gene}_expression_activity.png`：高相关性基因的表达量和基因活性 UMAP 图（如`PCDH9_expression_activity.png`、`NKAIN2_expression_activity.png`）
3. **分析结果表格**
    - 文件路径：`/data/project/r2py/young2-result/gene_atac_correlation.csv`
    - 格式：CSV 文件
    - 说明：RNA 表达与 ATAC 基因活性的相关性分析结果，包含基因名称、相关性系数和 p 值
4. **多组学数据结果**
    - `rna_atac_combined.h5mu`：完整的多组学分析结果，包含 RNA、ATAC 和基因活性三个模态
    - `rna_combined.h5ad`：处理后的 RNA 数据，包含联合分析结果（PCA、UMAP、聚类等）
    - `atac_combined.h5ad`：处理后的 ATAC 数据，包含联合分析结果

## **3.处理步骤**

1. **数据对齐**
    - 输入：RNA-seq 和 ATAC-seq 原始数据
    - 处理：保留两种测序中共同的细胞 ID
    - 输出：对齐后的多组学对象
    
    ```python
    # 细胞ID对齐
    common_cells = rna_adata.obs_names.intersection(atac_adata.obs_names)
    common_cells = list(common_cells)
    
    # RNA数据切片
    rna_adata = rna_adata[common_cells, :]
    
    # ATAC数据切片
    atac_adata = atac_adata.subset(
        obs_indices=common_cells,
        inplace=False
    )
    
    # 创建多组学对象
    mdata = mu.MuData({"rna": rna_adata, "atac": atac_adata})
    ```
    
2. **基因活性计算**
    - 输入：ATAC-seq 数据和基因组注释文件
    - 处理：使用 `snap.pp.make_gene_matrix` 计算基因活性矩阵
    - 输出：作为新模态（`gene_activity`）添加到多组学对象
    
    ```python
    # 计算基因活性矩阵
    gene_activity_adata = snap.pp.make_gene_matrix(
        mdata['atac'],
        gene_anno="/data/project/r2py/dataset/gencode_v41_GRCh38.gff3.gz",
        use_x=True,
        inplace=False,
    )
    ```
    
3. **联合降维和聚类**
    - 输入：RNA 表达矩阵和基因活性矩阵
    - 处理：通过标准化（`sc.pp.scale`）和特征合并（`np.hstack`）消除量纲差异。PCA 降维、UMAP 可视化、不同分辨率的聚类
    - 输出：PCA 方差解释率图、联合聚类 UMAP 图
    
    ```python
    # 提取RNA表达矩阵和ATAC基因活性矩阵
    rna_expr = mdata['rna'].X
    gene_activity = mdata['gene_activity'].X
    
    # 标准化处理
    rna_norm = sc.pp.scale(rna_expr, copy=True)
    activity_norm = sc.pp.scale(gene_activity, copy=True)
        
    # 特征匹配：找到两个数据集中共同的基因
    common_genes = list(set(mdata['rna'].var_names).intersection(set(mdata['gene_activity'].var_names)))
    
    # 提取共同基因的表达和活性
    rna_common_idx = [mdata['rna'].var_names.get_loc(gene) for gene in common_genes]
    activity_common_idx = [mdata['gene_activity'].var_names.get_loc(gene) for gene in common_genes]
        
    rna_common = rna_norm[:, rna_common_idx]
    activity_common = activity_norm[:, activity_common_idx]
        
    # 合并特征并进行PCA降维
    combined_features = np.hstack([rna_common, activity_common])
    pca = PCA(n_components=50)
    
    # 基于联合PCA结果进行UMAP和聚类
    sc.pp.neighbors(mdata['rna'], use_rep='X_joint_pca', n_neighbors=15)
    sc.tl.umap(mdata['rna'], min_dist=0.3)
        
    # 尝试不同分辨率的聚类
    resolutions = [0.3, 0.5, 0.8]
    for res in resolutions:
        key = f'joint_leiden_r{res}'
        sc.tl.leiden(mdata['rna'], resolution=res, key_added=key)
    ```
    
4. **关联分析**
    - 输入：RNA 表达和基因活性数据
    - 处理：计算高变基因的 RNA 表达与基因活性的 Spearman 相关性（Spearman 系数）
    - 输出：相关性表格和基因特异性的 UMAP 可视化
    
    ```python
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
    ```
    
5. **结果保存**
    - 输入：处理后的多组学对象
    - 处理：为每个模态的 `var_names` 添加唯一前缀，避免冲突
    - 输出：保存为`.h5mu` 和`.h5ad` 文件

## 4.输出结果

![joint_clustering.png](joint_clustering.png)

[analysis.log](analysis.log)

# 四、系统发育树（原R代码）

## 1.设置 redeemR 对象

### 1.1`redeemR.read`

**输入文件**：

- `QualifiedTotalCts` 存储每个细胞每个位置的四种线粒体 DNA 覆盖信息。
- `RawGenotypes.*.StrandBalance` 存储基因型信息和一致性水平等。

**输出文件**：

`VariantsGTSummary.RDS` 信息包括细胞条形码（Cell）、位置参考等位基因 / 变异等位基因（Variants）、突变等位基因数量（Freq）、总捕获等位基因数量（depth）、突变类型（type）、突变背景（Context）、异质性（hetero）

|  | **Var1** | **Cell** | **Variants** | **Freq** | **depth** | **Type** | **Context** | **hetero** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  | **<fct>** | **<chr>** | **<chr>** | **<int>** | **<dbl>** | **<chr>** | **<chr>** | **<dbl>** |
| **1** | AAACAAGCACCTCGAC_10626_T_G | AAACAAGCACCTCGAC | 10626_T_G | 8 | 26 | T_G | CTC | 0.30769231 |
| **2** | AAACAAGCACCTCGAC_11164_A_T | AAACAAGCACCTCGAC | 11164_A_T | 1 | 35 | A_T | CTA | 0.02857143 |
| **3** | AAACAAGCACCTCGAC_11545_T_A | AAACAAGCACCTCGAC | 11545_T_A | 1 | 37 | T_A | TTG | 0.02702703 |
- **源代码**
    
    ```r
    #' Function to generate GTS summary
    #'
    #' This function allows you to summarize the meta data for each genotyped variant
    #' @param RawGenotypes Well-named "RawGenotypes.Sensitive.StrandBalance" file in function redeemR.read or CW_mgatk.read
    #' @param filterN Boolean variable, if true filter out the variant with "N"
    #' @return Genotypes.summary a dataframe that summarize several metrics for each genotype
    #' @examples Usually used inside of function CW_mgatk.read
    #' @export
    #' @import dplyr
    GTSummary<-function(RawGenotypes,filterN=T){ ## At this moment, the context with N is probably prone to error due to mapping, in the future should work on realignment
    # 构建覆盖深度字典（Depthdic）
    data(ContextsDic)
    Depth<-unique(RawGenotypes[,c("Cell","Pos","Depth")])
    Depthdic<-Depth$Depth
    names(Depthdic)<-paste(Depth$Cell, Depth$Pos,sep="")
    # 统计突变出现次数，Cell和Variants
    Genotypes.summary<-table(paste(RawGenotypes$Cell,RawGenotypes$Variants,sep="_")) %>% as.data.frame()
    # 拆分字符串提取信息，从Cell和Variants的组合字符串（Var1）中拆分出Cell、Variants、cellPos
    Genotypes.summary$Cell<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){x[1]})
    Genotypes.summary$Variants<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[2:4],collapse="_")})
    Genotypes.summary$cellPos<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[1:2],collapse="")})
    # 补充关键指标
    Genotypes.summary$depth<-Depthdic[Genotypes.summary$cellPos]
    Genotypes.summary<-Genotypes.summary[,c("Var1","Cell","Variants","Freq","depth")]
    Genotypes.summary$Type<-strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){paste(x[2],x[3],sep="_")})
    Genotypes.summary$Context<-ContextsDic[strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){x[1]})]
    # 过滤低质量数据
    if(filterN){
        Genotypes.summary<-subset(Genotypes.summary,!grepl("N",Genotypes.summary$Context))
    }
    return(Genotypes.summary)
    }
    ```
    
- **python 代码**
    
    ```python
    def gt_summary(raw_genotypes, contexts_dic, filter_n=True):
        """
        生成每个基因型变体的元数据摘要，等效于R中的GTSummary函数
        
        参数:
            raw_genotypes: DataFrame，原始基因型数据（类似RawGenotypes.Sensitive.StrandBalance），需包含"Cell", "Pos", "Depth", "Variants"列
            contexts_dic: dict，位置-突变上下文字典（键为位置字符串，值为上下文序列）
            filter_n: bool，是否过滤Context中包含"N"的行
        
        返回:
            genotypes_summary: DataFrame，汇总后的基因型元数据
        """
        # 1. 构建覆盖深度字典（Depthdic）：键为Cell+Pos（无分隔符），值为Depth
        # 提取Cell、Pos、Depth的唯一组合
        depth_unique = raw_genotypes[["Cell", "Pos", "Depth"]].drop_duplicates()
        # 生成键：Cell和Pos拼接（无分隔符），值为Depth
        depth_dict = {f"{row['Cell']}{row['Pos']}": row['Depth'] for _, row in depth_unique.iterrows()}
        
        # 2. 统计每个Cell-Variants组合的出现次数（Freq）
        # 按Cell和Variants分组计数
        count_df = raw_genotypes.groupby(["Cell", "Variants"]).size().reset_index(name="Freq")
        # 生成Var1：Cell_Variants（用下划线连接）
        count_df["Var1"] = count_df["Cell"] + "_" + count_df["Variants"]
        
        # 3. 提取Pos（从Variants中拆分，Variants格式为"Pos_ref_alt"）
        count_df["Pos"] = count_df["Variants"].apply(lambda x: x.split("_")[0])
        
        # 4. 生成cellPos（Cell和Pos拼接，无分隔符），用于查询深度
        count_df["cellPos"] = count_df["Cell"] + count_df["Pos"]
        
        # 5. 补充depth（从depth_dict中查询）
        count_df["depth"] = count_df["cellPos"].map(depth_dict)
        
        # 6. 提取Type（突变类型，格式为"ref_alt"）
        count_df["Type"] = count_df["Variants"].apply(lambda x: "_".join(x.split("_")[1:3]))
        
        # 7. 补充Context（从contexts_dic中查询）
        count_df["Context"] = count_df["Pos"].map(contexts_dic)
        
        # 8. 过滤Context中包含"N"的行（若filter_n=True）
        if filter_n:
            # 排除Context含"N"或为空的行
            count_df = count_df[~count_df["Context"].str.contains("N", na=True)]
        
        # 9. 选择目标列并返回（与R结果列顺序一致）
        genotypes_summary = count_df[["Var1", "Cell", "Variants", "Freq", "depth", "Type", "Context"]]
        
        return genotypes_summary
    ```
    

### 1.2**`Create_redeemR`**

`redeemR` 对象充当一个容器，包含 ReDeeM 数据集的数据（如突变基因型）和分析结果（如克隆距离或系统发育结果）。

- `GTsummary.filtered` ：经过筛选的基因型信息，其中每一行代表给定细胞中的一个突变。

|  | **Var1** | **Cell** | **Variants** | **Freq** | **depth** | **Type** | **Context** | **hetero** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  | **<fct>** | **<chr>** | **<chr>** | **<int>** | **<dbl>** | **<chr>** | **<chr>** | **<dbl>** |
| **1** | AAACAAGCACCTCGAC_10626_T_G | AAACAAGCACCTCGAC | 10626_T_G | 8 | 26 | T_G | CTC | 0.30769231 |
| **2** | AAACAAGCACCTCGAC_11164_A_T | AAACAAGCACCTCGAC | 11164_A_T | 1 | 35 | A_T | CTA | 0.02857143 |
- `V.fitered` ：经过筛选的基因型特征信息，其中每一行代表所有细胞中的一个突变

|  | **Variants** | **CellN** | **PositiveMean** | **maxcts** | **CellNPCT** | **TotalVcount** | **TotalCov** | **totalVAF** | **CV** | **HomoTag** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  | **<chr>** | **<int>** | **<dbl>** | **<int>** | **<dbl>** | **<int>** | **<dbl>** | **<dbl>** | **<dbl>** | **<chr>** |
| **6** | 1000_T_C | 9 | 0.0328861 | 8 | 0.0014153169 | 22 | 451310 | 4.874698e-05 | 0.008177771 | Hetero |
| **11** | 10001_T_C | 3 | 0.0705272 | 12 | 0.0004717723 | 14 | 385839 | 3.628456e-05 | 0.122017337 | Hetero |
- **参数**
    - **Variants**突变的唯一标识，通常格式为 “位置_参考碱基_突变碱基”（如`1000_T_C`），表示 mtDNA 上第 1000 位的碱基由 T 突变为 C。
    - **CellN**该突变被检测到的细胞数量（Number of cells with the variant）。例如`CellN=9`表示该突变在 9 个细胞中出现过。
    - **PositiveMean**该突变在所有阳性细胞（即检测到该突变的细胞）中的平均异质性水平（mean heteroplasmy）。异质性（heteroplasmy）是指突变等位基因占该位置总等位基因的比例（即`hetero`列的值），`PositiveMean`反映该突变在不同细胞中的平均频率。
    - **maxcts**该突变在所有细胞中的最大覆盖深度（maximum depth across cells）。覆盖深度（depth）是指测序中覆盖该 mtDNA 位置的总 reads 数，`maxcts`代表该突变在单个细胞中被覆盖的最高次数。
    - **CellNPCT**携带该突变的细胞数占 “合格细胞总数” 的百分比（Percentage of qualified cells with the variant），计算公式为 `(CellN / QualifiedCellN) × 100`。例如`CellNPCT=0.0014`表示该突变仅在 0.14% 的合格细胞中出现。
    - **TotalVcount**该突变在所有细胞中的总突变等位基因数（total variant allele counts across cells），即所有细胞中`Freq`列的总和（`Freq`是单个细胞中该突变的等位基因数）。
    - **TotalCov**该突变在所有细胞中的总覆盖深度（total coverage across cells），即所有细胞中`depth`列的总和，反映该突变在所有细胞中的整体测序覆盖量。
    - **totalVAF**该突变的总变异等位基因频率（total Variant Allele Frequency），计算公式为 `TotalVcount / TotalCov`，反映该突变在所有细胞中的整体频率（而非单个细胞的异质性）。
    - **CV**变异系数（Coefficient of Variation），计算公式为 “该突变在不同细胞中异质性（hetero）的标准差 / PositiveMean”。CV 反映该突变的异质性在不同细胞中的波动程度：CV 越小，说明突变频率在细胞间越稳定；CV 越大，说明波动越剧烈。
    - **HomoTag**用于标记该突变是 “纯质性（Homoplasmy）” 还是 “异质性（Heteroplasmy）”：
        - 纯质性突变（`Homo`）：通常在大多数细胞中出现，且异质性接近 100%（即几乎所有该位置的等位基因都是突变型），可能是种系突变或早期克隆突变。
        - 异质性突变（`Hetero`）：仅在部分细胞中出现，或异质性较低（如 < 50%），通常是体细胞突变或后期克隆突变。
- **源代码**
    
    ```r
    #' Create_redeemR
    #'
    #' This function is to create redeemR with basic information
    #' @param VariantsGTSummary simply put GTSummary (Generated by redeemR.read) 
    #' @param qualifiedCellCut The minimum median mitochondrial coverage for a qualified cell, default is 10
    #' @param OnlyHetero If only consider the heteroplasmy variants, default is T
    #' @param VAFcut only use variants with VAF smaller than VAFcut. Default is 1.  We can use smaller value to constrain into only using rare variants
    #' @param Cellcut only use variants with at least cellcut cells carry
    #' @param maxctscut only use variants with at least in one cell with at leaset maxctscut variant fragments
    #' @return redeemR class
    Create_redeemR<-function(VariantsGTSummary=VariantsGTSummary,qualifiedCellCut=10,VAFcut=1,Cellcut=2,maxctscut=2){
    # 提取边缘修剪参数（edge_trim）,用于后续的序列边缘修剪
     if ("edge_trim" %in% names(attributes(VariantsGTSummary))){
            edge_trim <- as.numeric(attr(VariantsGTSummary,"edge_trim"))
        }else{
            edge_trim <- 0
        }
    # 筛选合格细胞（CellMeta）
    CellMeta<-subset(attr(VariantsGTSummary,"depth")[["Cell.MeanCov"]],meanCov>=qualifiedCellCut)
    names(CellMeta)[1]<-"Cell"
    # 过滤变异（Vfilter_v4），根据参数（最小细胞数Cellcut、单个细胞最大变异片段数maxctscut、合格细胞阈值qualifiedCellCut）过滤变异
    VariantsGTSummary.feature<-Vfilter_v4(VariantsGTSummary,Min_Cells = Cellcut, Max_Count_perCell = maxctscut, QualifyCellCut = qualifiedCellCut)
    # 生成过滤后的基因型数据（GTsummary.filtered），仅保留属于合格细胞且在过滤后变异列表中的基因型记录
    GTsummary.filtered<-subset(VariantsGTSummary,Variants %in% VariantsGTSummary.feature$Variants & Cell %in% CellMeta$Cell)
    # 构建redeemR对象
    ob<-new("redeemR")
    ob@GTsummary.filtered<-GTsummary.filtered
    ob@CellMeta<-CellMeta
    ob@V.fitered=VariantsGTSummary.feature
    ob@HomoVariants<-attr(VariantsGTSummary.feature,"HomoVariants")
    ob@UniqueV<-VariantsGTSummary.feature$Variants
    ob@DepthSummary<-attr(VariantsGTSummary,"depth")
    ob@para<-c(Threhold=attr(VariantsGTSummary,"thr"),qualifiedCellCut=qualifiedCellCut,VAFcut=VAFcut,Cellcut=Cellcut,maxctscut=maxctscut,edge_trim=edge_trim)
    ob@attr<-list(Filter.Cell=attr(VariantsGTSummary.feature,"Filter.Cell"),Filter.V=attr(VariantsGTSummary.feature,"Filter.V"),path=attr(VariantsGTSummary,"path"))
    return(ob)
    }
    ```
    

## 2.质量控制

- 线粒体 DNA 覆盖度（mtDNA coverage in each cell）
- 线粒体 DNA 突变特征（mtDNA mutation signature）
- 线粒体 DNA 突变丰度和频率（mtDNA mutation abundance and frequency）

## 3.线粒体 DNA 突变共识基准

`Show_Consensus`函数用于评估线粒体 DNA（mtDNA）突变的**一致性水平**，通过统计指标和可视化展示突变检测的可靠性。主要功能包括：

1. 读取对应严格度的原始基因型数据（`RawGenotypes.*.StrandBalance`）；
2. 筛选代表性变异（异质性、在多个细胞中出现）并抽样；
3. 计算关键指标（如一致性得分、链比例）；
4. 生成四幅可视化图表，展示 UMI 家族大小、一致性得分分布、测序重叠比例等；
5. 输出统计信息（分位数、重叠检测百分比）。

## 4.**细胞 - 突变关联矩阵**

`Make_matrix(Example_redeemR, onlyhetero=T)` ：将过滤后的 mtDNA 突变数据转换为**结构化的矩阵格式**，建立 “细胞” 与 “突变” 之间的关联关系。

- 参数`onlyhetero=T`表示仅保留**异质性突变**（Heteroplasmic mutations），因为纯质性突变（Homoplasmic mutations）在所有细胞中高度一致（如种系突变），无法用于区分不同克隆群（克隆关系依赖于细胞间突变的差异）。
- **生成稀疏矩阵：**两个存储在`redeemR`对象中的稀疏矩阵（`dgCMatrix`格式）：
    
    **（1）`@Cts.Mtx`：计数矩阵（Count Matrix）**
    
    - **含义**：矩阵的行是细胞（Cell），列是突变（Variants），矩阵中的值表示**某个细胞中该突变的 “突变等位基因数量”**（即`Freq`列的值，对应单个细胞中检测到的突变型碱基数量）。
    - **示例解读**：其中的`.`表示 “0”（稀疏矩阵的压缩表示方式），说明这 3 个细胞在这 3 个突变上均未检测到突变等位基因。
        
        ```
        3 x 3 sparse Matrix of class "dgCMatrix"  
                         Variants1000TC Variants10001TC Variants10002AC  
        AAACAAGCACCTCGAC              .               .               .  
        AAACCAGGTAACGAGC              .               .               .  
        AAACCAGGTGCAACGG              .               .               .
        ```
        
    
    **（2）`@Cts.Mtx.bi`：二值化矩阵（Binarized Matrix）**
    
    - **含义**：同样以细胞为行、突变为列，但矩阵中的值被二值化：
        - `1`：表示该细胞中检测到了该突变（无论突变等位基因数量多少，只要`Freq ≥ 1`）；
        - `0`：表示该细胞中未检测到该突变。
    - **作用**：简化分析，仅关注 “是否携带突变”，忽略突变数量的差异，适用于基于存在 / 缺失模式的克隆分型。
    - **为什么用 “稀疏矩阵”（`dgCMatrix`）？**
        
        单细胞 mtDNA 数据中，**每个细胞携带的突变数量远少于总突变数**（例如，总突变数可能有数千个，但单个细胞通常只携带几十个）。
        
        - 稀疏矩阵仅存储非零值，大幅减少内存占用和计算量（避免存储大量无意义的 0）；
        - `dgCMatrix`是 R 中`Matrix`包定义的高效稀疏矩阵格式，支持快速的矩阵运算（如后续的距离计算）。
- **源代码**
    
    ```r
    #' Make_matrix
    #' This will make the matixies of Cell VS mitochondrial variants and return redeemR
    #' Results stored in Cts.Mtx and Cts.Mtx.bi
    #' @param object redeemR class
    #' @param onlyhetero Only use heteroplasmic mutations
    #' @return redeemR class
    #' @export
    setMethod(f="Make_matrix",  # 定义方法名称
              signature="redeemR",  # 限定该方法仅适用于"redeemR"类对象
              definition=function(object, onlyhetero=T){  # 方法参数：redeemR对象、是否仅用异质性突变
                    require(dplyr)  # 加载数据处理包
                    require(Matrix.utils)  # 加载矩阵重塑包
                    
                    # 步骤1：筛选异质性突变（若参数为TRUE）
                    if(onlyhetero){
                        # 从过滤后的基因型数据中排除纯质性变异
                        GTsummary.filtered <- subset(object@GTsummary.filtered, !Variants %in% object@HomoVariants)
                        message("Only heteroplasmic mutations are used")  # 打印提示信息
                    }
                    
                    # 步骤2：将长格式数据重塑为细胞-突变计数矩阵
                    Cts.Mtx <- dMcast(GTsummary.filtered, Cell~Variants, value.var = "Freq")
                    
                    # 步骤3：统一矩阵列名格式
                    colnames(Cts.Mtx) <- strsplit(as.character(colnames(Cts.Mtx)), "_") %>% 
                        sapply(., function(x){paste(x[1], x[2], x[3], sep="")})
                    
                    # 步骤4：生成二值化矩阵
                    Cts.Mtx.bi <- Cts.Mtx  # 复制计数矩阵
                    Cts.Mtx.bi[Cts.Mtx.bi >= 1] <- 1  # 二值化：将"≥1"的数值转为1（表示存在突变）
                    
                    # 步骤5：将矩阵存储回redeemR对象
                    object@Cts.Mtx.bi <- Cts.Mtx.bi  # 存储二值化矩阵
                    object@Cts.Mtx <- Cts.Mtx  # 存储计数矩阵
                    message("@Cts.Mtx and @Cts.Mtx.bi are added")  # 打印提示信息
                    
                    return(object)  # 返回更新后的对象
    })
    
    ```
    

## 5.克隆距离与聚类

基于线粒体 DNA（mtDNA）突变数据推断细胞间克隆关系的完整流程，核心是通过**计算细胞间的突变距离**并进行**聚类 / 网络分析**，从而揭示细胞群体的克隆结构。

**距离计算→降维聚类→网络构建**

### **5.1核心思路：基于 mtDNA 突变的细胞关系推断**

mtDNA 突变可作为 “克隆标签”—— 同一克隆的细胞会共享相同的突变模式，而不同克隆的细胞突变模式差异较大。通过计算细胞间的 “突变距离”（差异程度），可以推断它们的克隆关系：距离越近，越可能属于同一克隆。

### **5.2距离计算的选择：二值化矩阵优先**

- 原因：`@Cts.Mtx`（计数矩阵）受测序深度偏差影响大（如高深度细胞可能检测到更多突变 allele），而`@Cts.Mtx.bi`（二值化矩阵）仅记录 “是否携带突变”，更稳定。
- 提供的距离方法：包括加权 Jaccard（w_jaccard）、LSI（潜在语义索引）、Jaccard、Dice 等，还支持 Simpson、Hamming 等其他二进制距离。

### **5.3LSI 降维与聚类（SeuratLSIClustering）**

这是基于 LSI 的聚类方法，流程如下：

1. **数据输入**：使用二值化矩阵`@Cts.Mtx.bi`，排除干扰性突变（如`rmvariants`参数指定的常见多态性）。
2. **标准化与降维**：
    - 用 TF-IDF（词频 - 逆文档频率）标准化数据，突出 “在部分细胞中高频出现的突变”（这些突变更可能是克隆标记）。
    - 用 SVD（奇异值分解）进行 LSI 降维，将高维突变数据压缩为 50 个核心维度（捕捉主要变异模式）。
3. **聚类与可视化**：
    - 用 UMAP/TSNE 将 LSI 结果降维到 2D 空间，便于可视化。
    - 用 Louvain 算法聚类（指定分辨率`res`），得到细胞亚群（潜在克隆群）。
4. **结果**：生成 Seurat 对象存储在`@Seurat`中，可通过 UMAP 图展示聚类结果（不同颜色代表不同克隆群）。
- **源代码**
    
    ```r
    #' SeuratLSIClustering
    #' This will use the mito variants for Seurat clustering (LSI based)
    #' @param  redeemR class
    #' @param binary  Default is tree, to make use of the binary matrix
    #' @param res     Default os 0.3, the resolution of the clustering
    #' @return redeemR class
    #' @export
    setMethod(f="SeuratLSIClustering",
              signature="redeemR",
              definition=function(object,binary=T,res=0.6,lsidim=2:50,rmvariants=c("Variants310TC","Variants3109TC","Variants5764CT")){
              require(Signac)
              require(Seurat)
              # 数据准备与过滤
              if(binary){
                  if (packageVersion("Seurat")>"4.0.0"){
                      print("Seurat5 is on, convert data structure to v3")
                      options(Seurat.object.assay.version = 'v3')
                  }
                  Cts.Mtx.bi<-as.matrix(object@Cts.Mtx.bi)
                  # 移除rmvariants中指定的变异
                  Cts.Mtx.bi<-Cts.Mtx.bi[,!colnames(Cts.Mtx.bi) %in% rmvariants]
                  # 过滤无突变的细胞
                  Cts.Mtx.bi<-Cts.Mtx.bi[rowSums(Cts.Mtx.bi)>0,]
                  # 创建 Seurat 对象
                  Cell_Variant.seurat<-CreateSeuratObject(counts = as(t(as.matrix(Cts.Mtx.bi)),"CsparseMatrix"), assay = "redeemR")
              }else{
                  Cts.Mtx<-as.matrix(object@Cts.Mtx)
                  Cts.Mtx<-Cts.Mtx[,!colnames(object@Cts.Mtx) %in% rmvariants]
                  Cts.Mtx<-Cts.Mtx[rowSums(Cts.Mtx)>0,]
                  Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(Cts.Mtx)), assay = "redeemR")
              }
              
              # 特征处理与标准化
              # 设置可变特征，将所有变异设为 “可变特征”
              VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
              # TF-IDF 转换，对数据进行标准化，突出在部分细胞中高频出现的突变
              Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
              
              # LSI 降维
              # 筛选顶级特征，保留至少在 2 个细胞中出现的突变
              Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 2)
              # 奇异值分解（SVD），将高维的细胞 - 突变矩阵压缩为 50 个潜在维度
              Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
              
              # 聚类与可视化
              # RunUMAP和RunTSNE基于 LSI 结果进行非线性降维，将细胞映射到 2D 空间
              Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = lsidim)
              Cell_Variant.seurat <- RunTSNE(Cell_Variant.seurat, reduction = "lsi", dims = lsidim,check_duplicates = FALSE)
              # FindNeighbors基于 LSI 维度计算细胞间距离，FindClusters以指定分辨率（res）进行聚类，得到细胞亚群
              Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = lsidim)
              Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = res)
              object@Seurat<-Cell_Variant.seurat
              return(object)
    })
    ```
    
- python代码
    - 使用 PCA 替代 LSI（在 scanpy 中没有直接对应的 LSI 实现，但 PCA 可以实现类似的降维效果）
    - 使用 Leiden 聚类算法替代 Seurat 的 Louvain 算法

### **5.4加权距离：校正平行进化的影响**

- **平行进化问题**：不同克隆的细胞可能独立出现相同突变（非遗传导致），会干扰距离计算。
- **解决方案**：通过 “权重” 校正 —— 基于群体中突变的复发率（`CellPCT`数据），给 “高频复发突变”（更可能平行进化）更低权重，给 “罕见突变”（更可能是克隆特异性标记）更高权重。
- **操作**：
    - 准备权重数据`V.weight`（`weight=1-突变复发率`）。
    - 用`AddDist`函数计算并存储多种加权距离（如 w_jaccard、w_cosine）和非加权距离（如 Jaccard、Dice），结果存在`@DistObjects`中。
- **源代码**
    
    ```r
    #' AddDist
    #' This add Jaccard, Dice, Jaccard3W distance and stored in DistObjects
    #' @param object redeemR class
    #' @param jaccard  default=T
    #' @param dice    default=T
    #' @param jaccard3w  default=T
    #' @param w_jaccard   default=T
    #' @param w_cosine default=T
    #' @param weight A two column dataframe, "Variant"(The variant name should match cell-variant matrix column, e.g, Variants310TC), "weight" (numeric)
    #' @param NN To replace NA, which means a variant shown in the object is not shown in the weight vector, with a number, default is 1 for jaccard system. 
    #' @param LSIdist default=T
    #' @param dim the dimensions to use to calculate LSI distance default is 2:50
    #' @return redeemR class
    #' @export
    setMethod(f="AddDist",
              signature="redeemR",
              definition=function(object,jaccard=T,dice=T,jaccard3w=T,w_jaccard=T,w_cosine=T,weightDF=NULL,NN=1,LSIdist=T,dim=2:50){
              d.Jaccard<-NA
              d.Dice<-NA    
              d.3WJaccard<-NA
              d.w_jaccard<-NA
              d.w_cosine<-NA
              
              # 权重处理
              if(length(weightDF)!=0){
                weight<-data.frame(Variants=colnames(object@Cts.Mtx.bi)) %>% merge(.,weightDF,by="Variants",all.x = T,sort = F) %>% .$weight
              }  
              if(length(which(is.na(weight)))!=0){
                weight[is.na(weight)]<-NN
                print("Some variant i weight is not found in cell-variant matrix, use 1")
              }
              if(length(weight)!=ncol(object@Cts.Mtx.bi)){
                 stop("The length of weight does not match the variant numbers in the martix")
              }
              print("Weight vector matches well with the Cell-Variant matrix, continue...")
              
              # 距离计算
              if(jaccard){
                  d.Jaccard<-BinaryDist(object@Cts.Mtx.bi,method="Jaccard")
                  message("jaccard distances added")
              }    
              if(dice){
                   d.Dice<-BinaryDist(object@Cts.Mtx.bi,method="Dice")
                   message("dice distances added")
              }
              if(jaccard3w){
                  d.3WJaccard<-BinaryDist(object@Cts.Mtx.bi,method="3WJaccard")
                  message("3wjaccard distances added")
              }
              if(w_jaccard){
                  if(length(weightDF)==0){
                      stop("Please input the weight, otherwise turn off the w_jaccard")
                      
                  }
                  d.w_jaccard<-quick_w_jaccard(object@Cts.Mtx.bi,w=weight)
                  message("weighted jaccard distances added")
              }
              if(w_cosine){
                  if(length(weightDF)==0){
                      stop("Please input the weight, otherwise turn off the w_cosine")
                  }
                  d.w_cosine<-quick_w_cosine(object@Cts.Mtx.bi,w=weight)
                  message("weighted cosine distances added")
              }
              if(LSIdist){
                  d.lsi<-dist(object@Seurat@reductions$lsi@cell.embeddings[,dim])
                  message("LSI distances added")
              }
              object@DistObjects<-new("DistObjects",jaccard=d.Jaccard, Dice=d.Dice,jaccard3W=d.3WJaccard,w_jaccard=d.w_jaccard,w_cosine=d.w_cosine,LSIdist=d.lsi)
              return(object)
              })
    ```
    

### **5.5细胞 - 细胞网络构建：基于距离的互近邻（MNN）网络**

通过距离矩阵构建细胞间的相似性网络，同一克隆的细胞在网络中连接更紧密：

1. **FromDist2Graph**：直接从距离矩阵（如`@DistObjects@w_jaccard`）生成 igraph 对象（网络），节点是细胞，边代表细胞间的相似性（距离越近，连接越可能存在）。
2. **KNN 邻域方法**：
    - `MakeNN`：为每个细胞找最近的 k 个邻居（如 15 个），输出邻居索引和距离。
    - `NN2M`：将 KNN 结果转换为邻接矩阵。
    - `graph_from_adjacency_matrix`：将邻接矩阵转换为 igraph 网络。
- **源代码**
    
    ```r
    #' FromDist2Graph 
    #' From distance object or matrix to graph, default is to return igraph object
    #' This function was developed based on 
    #' @param d the distance matrix,  this can be either dist or a matrix
    #' @param k.param K default is 30
    #' @param return_igraph Wheather return igraph, default is T which return igraph. Otherwise, return adjacent matrix
    #' @return igraph or adjacent matrix
    #' @export
    #' @import Matrix
    #' @importFrom  igraph get.adjacency graph.edgelist graph_from_adjacency_matrix
    FromDist2Graph<-function(d,k.param=30,return_igraph=T){
    # 距离矩阵格式统一
    if(!is.matrix(d)){
      d<-as.matrix(d)  # 若输入是dist对象，转换为矩阵
    }
    # 初始化 KNN 存储矩阵
    n.cells <- dim(d)[1]  # 细胞数量（矩阵行数）
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)  # 存储每个细胞的K个最近邻索引
    knd.mat <- knn.mat  # 存储每个细胞的K个最近邻的距离值
    # 计算每个细胞的 K 最近邻（KNN）
    for (i in 1:n.cells) { 
      knn.mat[i, ] <- order(d[i, ])[1:k.param]  # 对第i个细胞的距离排序，取前k.param个索引（最近的K个细胞）
      knd.mat[i, ] <- d[i, knn.mat[i, ]]  # 提取这K个最近邻对应的距离值
    } 
    nn.dist <- knn.mat[, 1:k.param]  # 简化为仅保留最近邻索引（变量名可能有误，实际存储的是索引）
    knn <- nn.dist  # 重命名为knn，明确存储的是每个细胞的K个最近邻索引
    # nn.dist<- knd.mat[, 1:k.param]
    # nn<-list(idx=nn.idx,dist=nn.dist)
    # Compute MNN  -- Borrowed from SCAVENGE https://github.com/sankaranlab/SCAVENGE
    # 筛选互近邻（MNN）
    knn2 <- list()  # 存储每个细胞的互近邻边
    length(knn2) <- nrow(knn)  # 列表长度等于细胞数量
    
    for (i in 1:nrow(knn)) {
      # 检查细胞i的K个近邻（排除自身，取knn[i, -1]）中，哪些细胞的近邻包含i
      xx <- apply(knn[knn[i, -1], ], 1, function(x) {
        any(x == i)  # 判断细胞i是否是该近邻的近邻
      })
      
      if (sum(xx) > 0) {
        # 若存在互近邻，保留这些双向连接
        temp_knn <- knn[i, c(TRUE, xx)]  # 保留自身和互近邻
        temp_el <- cbind(temp_knn[1], c(temp_knn[-1]))  # 构建边（i与互近邻的连接）
      } else {
        # 若没有互近邻，保留与最近的1个细胞的连接（避免孤立节点）
        temp_el <- knn[i, 1:2]
      }
      knn2[[i]] <- temp_el  # 存储当前细胞的互近邻边
    }
    # 构建网络结构
    el <- do.call(rbind.data.frame, knn2) %>% as.matrix  # 将所有边整合为边列表（每行是一条边：两个细胞的索引）
    adj <- igraph::get.adjacency(igraph::graph.edgelist(el))  # 从边列表生成邻接矩阵
    mutualknn <- 1 * ((adj + Matrix::t(adj)) > 0)  # 确保邻接矩阵是无向的（对称），值为1表示存在连接
    colnames(mutualknn) <- rownames(mutualknn) <- rownames(d)  # 用细胞名命名行和列
    # 返回结果
    if(return_igraph){
      g <- igraph::graph_from_adjacency_matrix(mutualknn, diag = F, mode = "undirected")  # 转换为igraph对象
      return(g)
    }else{
      return(mutualknn)  # 返回邻接矩阵
    }
    }
    ```
    

## 6.单细胞系统发育树分析

### **6.1构建系统发育树**

基于已计算的细胞间距离，使用邻接法（neighbor joining, NJ）构建系统发育树。

通过`Make_tree`函数实现，构建的树存储在对象`Example_redeemR`的`@TREE`属性中。

- **源代码**
    
    ```r
    #' Make_tree
    #' This will generate a basic phylogenetic tree
    #' @param object redeemR class
    #' @param d "jaccard" or "Dice" or "jaccard3W" or  "w_jaccard"  "w_cosine"  "LSIdist"
    #' @param algorithm the algorithm used to build the tree, choose from "nj" and "upgma"
    #' @return redeemR class
    #' @export
    setMethod(f="Make_tree",
              signature="redeemR",
              definition=function(object,d,algorithm,onlyreturntree=F){
              dist<-slot(object@DistObjects,d)
              # 选择算法构建树
              if(algorithm=="nj"){
              phylo<-nj(dist)
              }else if (algorithm=="upgma"){
              phylo<-upgma(dist)
              }
              # 转换树格式
              treedata<-as.treedata(phylo)
              # 创建 TREE 对象，实例化一个自定义的TREE类对象
              TREEobject<-new("TREE",phylo=phylo,treedata=treedata,records=paste(d,algorithm,sep="-"))
              if(onlyreturntree){
              return(TREEobject)
              }else{
              object@TREE<-TREEobject
              return(object)
              }
    })
    ```
    

### **6.2可视化系统发育树**

使用`ggtree`函数可视化构建的树。

### **6.3更新细胞元数据（CellMeta）**

为细胞元数据添加两个关键信息：线粒体突变的变异数（VN）和细胞在树上的位置（TreePos）：

- 从树的可视化结果中提取细胞位置信息（包括细胞 ID、节点编号、y 坐标），命名为`CellTreePos`；
- 计算每个细胞的线粒体突变数（`VN`），通过对`Example_redeemR@Cts.Mtx.bi`矩阵按行求和得到，并与`CellTreePos`合并为`VN.summary`；
- 将`VN.summary`中的信息（VN、node、TreePos）合并到`Example_redeemR@CellMeta`中（若`TreePos`不存在于元数据中）。

# 五、系统发育树（python代码）

[young2_buildTree.zip](young2_buildTree.zip)

## 1.目录结构

```python
├── young2_buildTree
│   └── src
│       ├── mito_consensus.py  # 用于处理 mtDNA 突变数据的工具函数集
│       ├── dist_cal.py  #  围绕 mtDNA 突变数据的距离计算与网络构建的工具函数集
│       ├── main.py  # 实现从 mtDNA 突变数据的读取、预处理，到聚类分析、距离计算、网络构建和系统发育树推断的完整分析流程
│       ├── objects.py  # 用于 mtDNA 突变数据的聚类分析与系统发育树构建的工具类
│       ├── tree_between_clusters.py  # 基于 mtDNA 突变数据的聚类簇代表细胞选择与系统发育树构建的分析流程 
│       ├── visual.py  # 绘制圆形布局的系统发育树
│       ├── result
│       │   ├── cluster_representative_analysis  # 聚类簇代表细胞输出结果
│       │   │   ├── cluster_representative_circular_tree.png  # 聚类簇代表圆形布局细胞系统发育树
│       │   │   ├── cluster_representative_tree.png  # 聚类簇代表方形布局细胞系统发育树
│       │   │   └── representative_tree_result.pkl
│       │   ├── example_tree.png
│       │   ├── lsi_clustering_umap.png
│       │   ├── mnn_network.png
│       │   └── redeemr_with_dist.pkl
│       └── __pycache__
```

## 2.输出结果

全部细胞（因为构建树耗时太久没有发育树结果）

`mnn_network.png` ：

![mnn_network.png](mnn_network.png)

`lsi_clustering_umap.png` ：

![lsi_clustering_umap.png](lsi_clustering_umap.png)

聚类簇代表细胞（简化为79个细胞进行发育树构建）

`cluster_representative_circular_tree.png`

![cluster_representative_circular_tree.png](cluster_representative_circular_tree.png)

`cluster_representative_tree.png`

![cluster_representative_tree.png](cluster_representative_tree.png)
