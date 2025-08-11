import os
import snapatac2 as snap
import scanpy as sc
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt

# 辅助函数：验证并转换元数据列类型
def validate_and_convert_column(df, col_name):
    """验证DataFrame列的数据类型是否被AnnData支持，并在必要时进行转换"""
    dtype = df[col_name].dtype
    
    # 处理数值类型
    if pd.api.types.is_numeric_dtype(dtype):
        return df[col_name]
    
    # 处理字符串类型
    elif pd.api.types.is_string_dtype(dtype) or pd.api.types.is_object_dtype(dtype):
        # 检查是否可以转换为数值
        try:
            numeric_col = pd.to_numeric(df[col_name], errors='raise')
            print(f"列 {col_name}: 从字符串转换为数值类型")
            return numeric_col
        except (ValueError, TypeError):
            # 转换为分类类型
            print(f"列 {col_name}: 转换为分类类型")
            return df[col_name].astype('category')
    
    # 处理布尔类型
    elif pd.api.types.is_bool_dtype(dtype):
        return df[col_name]
    
    # 处理分类类型
    elif pd.api.types.is_categorical_dtype(dtype):
        return df[col_name]
    
    # 其他类型默认转换为字符串分类
    else:
        print(f"列 {col_name}: 未知类型 ({dtype})，转换为分类类型")
        return df[col_name].astype(str).astype('category')

# 设置中文字体，确保中文正常显示
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei"]

# 定义输出目录
output_dir = "/data/project/r2py/young2ATAC-result"
os.makedirs(output_dir, exist_ok=True)

# 读取数据
print("读取数据文件...")
try:
    # 读取稀疏矩阵并转置
    counts = scipy.io.mmread("/data/project/r2py/young2ATAC_need/atac_counts.mtx").tocsr().T
    print(f"读取的稀疏矩阵形状: {counts.shape}")
    
    # 读取peak注释和细胞元数据
    peaks = pd.read_csv("/data/project/r2py/young2ATAC_need/peak_annotations.csv", index_col=0)
    meta = pd.read_csv("/data/project/r2py/young2ATAC_need/cell_metadata.csv", index_col=0)
    
    print(f"Peak注释数量: {len(peaks)}")
    print(f"细胞元数据数量: {len(meta)}")
    
    # 确保细胞数量一致
    if counts.shape[0] != len(meta):
        print(f"警告: 矩阵中的细胞数({counts.shape[0]})与元数据中的细胞数({len(meta)})不一致!")
    
    # 确保peak数量一致
    if counts.shape[1] != len(peaks):
        print(f"警告: 矩阵中的peak数({counts.shape[1]})与peak注释数({len(peaks)})不一致!")
    
except Exception as e:
    print(f"读取数据时出错: {e}")
    exit(1)

# 创建AnnData对象
print("创建AnnData对象...")
try:
    # 指定HDF5文件路径
    h5ad_file = os.path.join(output_dir, "initial_data.h5ad")
    
    # 确保counts是scipy CSR矩阵
    if not isinstance(counts, scipy.sparse.csr_matrix):
        print("将counts矩阵转换为CSR格式...")
        counts = counts.tocsr()  # 强制转换为CSR
    
    # 这里counts的行是细胞，列是peak
    adata = snap.AnnData(filename=h5ad_file)
    adata.X = counts
    # 使用更安全的方式设置obs和var
    adata.obs_names = meta.index
    adata.var_names = peaks.index
    
    # 添加观察值(obs)，验证并转换每一列
    print("添加细胞元数据...")
    for col in meta.columns:
        # 检查是否有缺失值
        if meta[col].isna().any():
            print(f"警告: 列 {col} 包含缺失值，将填充为适当的默认值")
            if pd.api.types.is_numeric_dtype(meta[col]):
                meta[col] = meta[col].fillna(0)
            elif pd.api.types.is_string_dtype(meta[col]) or pd.api.types.is_object_dtype(meta[col]):
                meta[col] = meta[col].fillna('nan')
        
        # 验证并转换列类型
        converted_col = validate_and_convert_column(meta, col)
        adata.obs[col] = converted_col
    
    # 添加变量(var)，验证并转换每一列
    print("添加peak注释...")
    for col in peaks.columns:
        # 检查是否有缺失值
        if peaks[col].isna().any():
            print(f"警告: 列 {col} 包含缺失值，将填充为适当的默认值")
            if pd.api.types.is_numeric_dtype(peaks[col]):
                peaks[col] = peaks[col].fillna(0)
            elif pd.api.types.is_string_dtype(peaks[col]) or pd.api.types.is_object_dtype(peaks[col]):
                peaks[col] = peaks[col].fillna('nan')
        
        # 验证并转换列类型
        converted_col = validate_and_convert_column(peaks, col)
        adata.var[col] = converted_col
        
    print(adata)
    
    # 添加峰坐标信息（如果peaks DataFrame包含chr、start、end列）
    if 'chr' in peaks.columns and 'start' in peaks.columns and 'end' in peaks.columns:
        adata.var['chrom'] = peaks['chr']
        adata.var['start'] = peaks['start']
        adata.var['end'] = peaks['end']
    
    print(f"AnnData对象创建完成，包含 {adata.n_obs} 个细胞和 {adata.n_vars} 个特征")
    
except Exception as e:
    print(f"创建AnnData对象时出错: {e}")
    exit(1)

"""
step1：降维
"""
print("准备下游分析数据...")
try:
    # 直接基于 peak 数据选择高变特征
    snap.pp.select_features(adata, n_features=250000)  # 基于 peak 计数选择特征
    
except Exception as e:
    print(f"数据预处理时出错: {e}")
    
print("执行降维分析...")
try:
    # 谱嵌入降维（保留数据全局结构）
    snap.tl.spectral(adata)
    
    # UMAP降维
    snap.tl.umap(adata)
    
    print("降维完成，结果保存在 data.obsm['X_umap']")
    
    # 可视化降维结果
    snap.pl.umap(
        adata,
        color=['n_counts', 'tsse'],
        out_file=os.path.join(output_dir, "umap_initial.png"),
        height=500,
        width=1200,
    )
except Exception as e:
    print(f"降维时出错: {e}")

"""
step2：聚类分析
"""
print("执行聚类分析...")
try:
    # 构建K近邻图
    snap.pp.knn(adata)
    
    # Leiden算法聚类（基于图结构）
    snap.tl.leiden(adata)
    
    # 绘制UMAP聚类图
    snap.pl.umap(
        adata,
        color='leiden',
        out_file=os.path.join(output_dir, "umap_leiden.png"),
        height=500,
        width=600,
    )
    
    print(f"聚类完成，共识别出 {len(adata.obs['leiden'].unique())} 个细胞簇")
except Exception as e:
    print(f"聚类时出错: {e}")

