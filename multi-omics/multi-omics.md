# 一、`young2_RNA.py`

[young2_RNA.py](young2_ATAC.py)

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


