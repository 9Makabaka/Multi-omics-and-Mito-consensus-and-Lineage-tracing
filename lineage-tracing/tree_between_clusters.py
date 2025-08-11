import os
import pickle
import numpy as np
from matplotlib import pyplot as plt
from mito_consensus import create_redeemr, make_matrix, generate_cell_pct, read_rds_variants
from dist_cal import add_dist
from Bio.Phylo import draw
from visual import draw_circular_tree

plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]  # 支持中文的字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# --------------------------
# 1. 加载数据和聚类结果
# --------------------------
save_path_dist = "/data/project/r2py/young2_mitoConsensus/src/result/redeemr_with_dist.pkl"
with open(save_path_dist, 'rb') as f:
    redeemr = pickle.load(f)
print("redeemr对象加载完成")

# 检查聚类结果
if redeemr.seurat is None or 'leiden' not in redeemr.seurat.obs:
    raise ValueError("未找到聚类结果，请先运行LSI聚类分析")

# --------------------------
# 2. 为每个簇选择代表细胞（以“中心细胞”为例）
# --------------------------
def select_representative_cells(redeemr, cluster_col='leiden'):
    """选择簇内中心细胞作为代表（优化：增加距离矩阵检查）"""
    cluster_labels = redeemr.seurat.obs[cluster_col]
    all_cells = cluster_labels.index.tolist()
    cell_to_idx = {cell: i for i, cell in enumerate(all_cells)}
    
    # 检查距离矩阵是否存在
    if not hasattr(redeemr.dist_objects, 'w_jaccard'):
        raise AttributeError("redeemr.dist_objects中未找到w_jaccard距离矩阵，请先计算")
    dist_matrix = redeemr.dist_objects.w_jaccard
    
    representatives = {}
    for cluster_id in cluster_labels.unique():
        cluster_cells = cluster_labels[cluster_labels == cluster_id].index.tolist()
        if len(cluster_cells) < 1:
            print(f"警告：簇{cluster_id}无细胞，跳过")
            continue
        
        # 转换细胞ID为索引
        try:
            cluster_indices = [cell_to_idx[cell] for cell in cluster_cells]
        except KeyError as e:
            print(f"警告：簇{cluster_id}中存在未知细胞ID {e}，跳过该细胞")
            cluster_indices = [cell_to_idx[cell] for cell in cluster_cells if cell in cell_to_idx]
            if not cluster_indices:
                print(f"警告：簇{cluster_id}无有效细胞，跳过")
                continue
        
        # 计算平均距离（处理空矩阵情况）
        avg_distances = []
        for idx in cluster_indices:
            try:
                dists = dist_matrix[idx, cluster_indices]
                avg_dist = np.mean(dists)
                avg_distances.append((idx, avg_dist))
            except IndexError:
                print(f"警告：簇{cluster_id}中细胞索引{idx}超出距离矩阵范围，跳过")
                continue
        
        if not avg_distances:
            print(f"警告：簇{cluster_id}无法计算平均距离，随机选择一个细胞作为代表")
            min_idx = cluster_indices[0] if cluster_indices else None
        else:
            min_idx, _ = min(avg_distances, key=lambda x: x[1])
        
        if min_idx is not None:
            representative_cell = all_cells[min_idx]
            representatives[cluster_id] = representative_cell
            print(f"簇 {cluster_id} 的代表细胞：{representative_cell}（簇内细胞数：{len(cluster_cells)}）")
    
    return representatives

# 执行代表细胞选择
representatives = select_representative_cells(redeemr, cluster_col='leiden')
representative_cells = list(representatives.values())
if not representative_cells:
    raise ValueError("未成功选择任何代表细胞，请检查聚类结果")

# --------------------------
# 3. 基于代表细胞构建系统发育树（核心优化）
# --------------------------
def build_tree_with_representatives(redeemr, representative_cells, representatives, output_dir):
    """优化：增强ID匹配、增加调试信息、容错处理"""
    # 复用原始variants_gtsummary，避免重复读取导致差异
    rds_file_path = "/data/dataset/young2_HSC/Youn2.HSC.Consensus.final/VariantsGTSummary.RDS"
    
    # 步骤1：读取RDS文件（提取"Sensitive"严格度的数据）
    print(f"读取文件: {rds_file_path}，提取严格度：S（Sensitive）")
    variants_gtsummary = read_rds_variants(
        rds_path=rds_file_path,
        thr="S"
    )
    
    # 创建rep_redeemr（确保参数与原始一致）
    try:
        # 从原始对象提取参数（兼容不同存储方式）
        params = {
            'qualified_cell_cut': redeemr.qualified_cell_cut if hasattr(redeemr, 'qualified_cell_cut') else redeemr.para.get('qualified_cell_cut'),
            'vaf_cut': redeemr.vaf_cut if hasattr(redeemr, 'vaf_cut') else redeemr.para.get('vaf_cut'),
            'cell_cut': redeemr.cell_cut if hasattr(redeemr, 'cell_cut') else redeemr.para.get('cell_cut'),
            'maxcts_cut': redeemr.maxcts_cut if hasattr(redeemr, 'maxcts_cut') else redeemr.para.get('maxcts_cut'),
            'thr': redeemr.thr if hasattr(redeemr, 'thr') else redeemr.para.get('thr', 'S')
        }
        # 过滤None值（使用默认值）
        params = {k: v for k, v in params.items() if v is not None}
        
        rep_redeemr = create_redeemr(
            variants_gtsummary=variants_gtsummary,** params
        )
    except Exception as e:
        raise RuntimeError(f"创建rep_redeemr失败：{str(e)}")
    
    # 生成矩阵
    rep_redeemr = make_matrix(rep_redeemr, onlyhetero=True)
    
    # 检查cts_mtx有效性
    if rep_redeemr.cts_mtx is None or 'rows' not in rep_redeemr.cts_mtx:
        raise ValueError("rep_redeemr.cts_mtx初始化失败，请检查make_matrix函数")
    cts_rows = rep_redeemr.cts_mtx['rows']
    cts_matrix = rep_redeemr.cts_mtx['matrix']
    print(f"\nrep_redeemr.cts_mtx包含{len(cts_rows)}个细胞，矩阵形状：{cts_matrix.shape}")
    
    # 规范化ID匹配（处理可能的格式差异）
    def normalize_id(cell_id):
        """根据实际数据格式调整：去除特殊字符/统一大小写"""
        return cell_id.strip().upper()  # 示例：去除空格并转为大写
    
    # 构建规范化ID映射
    cts_id_map = {normalize_id(cell): cell for cell in cts_rows}
    valid_reps = []
    missing_reps = []
    for rep in representative_cells:
        norm_rep = normalize_id(rep)
        if norm_rep in cts_id_map:
            valid_reps.append(cts_id_map[norm_rep])
        else:
            missing_reps.append(rep)
    
    # 打印匹配结果
    print(f"\n代表细胞总数：{len(representative_cells)}")
    print(f"匹配成功的代表细胞：{len(valid_reps)}个")
    if missing_reps:
        print(f"匹配失败的代表细胞（前5个）：{missing_reps[:5]}")
    
    # 容错：至少保留3个有效细胞
    if len(valid_reps) < 3:
        if len(valid_reps) == 0:
            raise ValueError("无匹配的代表细胞，无法构建树")
        print(f"警告：有效代表细胞不足3个（{len(valid_reps)}个），可能影响树的可靠性")
    
    # 筛选矩阵（只保留有效代表细胞）
    try:
        valid_indices = [cts_rows.index(cell) for cell in valid_reps]
        # 更新rows和matrix
        rep_redeemr.cts_mtx['rows'] = [cts_rows[i] for i in valid_indices]
        rep_redeemr.cts_mtx['matrix'] = cts_matrix[valid_indices, :]
        print(f"筛选后矩阵形状：{rep_redeemr.cts_mtx['matrix'].shape}")
        
        if rep_redeemr.cts_mtx_bi is not None and 'matrix' in rep_redeemr.cts_mtx_bi:
            # 确保二值化矩阵的rows与原始矩阵一致（否则重新对齐）
            if 'rows' in rep_redeemr.cts_mtx_bi:
                rep_redeemr.cts_mtx_bi['rows'] = [rep_redeemr.cts_mtx_bi['rows'][i] for i in valid_indices]
            # 筛选二值化矩阵的行（仅保留代表细胞）
            rep_redeemr.cts_mtx_bi['matrix'] = rep_redeemr.cts_mtx_bi['matrix'][valid_indices, :]
            print(f"筛选后cts_mtx_bi形状：{rep_redeemr.cts_mtx_bi['matrix'].shape}")
        else:
            # 如果二值化矩阵不存在，尝试重新生成（避免后续计算出错）
            print("警告：cts_mtx_bi不存在，重新生成二值化矩阵...")
            rep_redeemr = make_matrix(rep_redeemr, onlyhetero=True)  # 基于筛选后的cts_mtx重新生成
            if rep_redeemr.cts_mtx_bi is not None:
                print(f"重新生成的cts_mtx_bi形状：{rep_redeemr.cts_mtx_bi['matrix'].shape}")
            else:
                raise ValueError("无法生成cts_mtx_bi，后续距离计算可能失败")
    except Exception as e:
        raise RuntimeError(f"筛选cts_mtx失败：{str(e)}")
    
    rep_redeemr.dist_objects = None  # 强制清空旧的距离矩阵
    print("已重置dist_objects，将基于筛选后的细胞重新计算距离")
    # 为筛选后的代表细胞重新运行LSI聚类（基于筛选后的矩阵）
    print("\n为筛选后的代表细胞运行LSI聚类...")
    try:
        # 获取筛选后矩阵的维度
        n_cells = len(rep_redeemr.cts_mtx['rows'])
        n_variants_original = rep_redeemr.cts_mtx['matrix'].shape[1]
        print(f"筛选后矩阵维度：细胞数={n_cells}, 原始变异数={n_variants_original}")
    
        # 检查是否满足min(细胞数, 变异数) > 50（因为函数固定n_comps=50）
        required_min_dim = 51  # 必须满足 50 < min(...)
        current_min_dim = min(n_cells, n_variants_original)
    
        # 如果变异数不足，调整变异筛选策略（保留更多变异）
        if current_min_dim < required_min_dim:
            # 检查是否是变异数不足导致
            if n_variants_original < required_min_dim and n_cells >= required_min_dim:
                print(f"变异数不足（{n_variants_original}），尝试保留更多变异...")
                # 重新生成矩阵时不过滤太多变异（仅过滤必要的）
                rep_redeemr = make_matrix(rep_redeemr, onlyhetero=True, min_cell_cut=1)  # 降低细胞数过滤阈值
                # 重新获取变异数
                n_variants_updated = rep_redeemr.cts_mtx['matrix'].shape[1]
                print(f"调整后变异数：{n_variants_updated}")
                current_min_dim = min(n_cells, n_variants_updated)
        
            # 如果细胞数不足，提示无法解决（需要更多代表细胞）
            if n_cells < required_min_dim:
                raise ValueError(f"细胞数不足（{n_cells}），至少需要{required_min_dim}个代表细胞才能运行LSI聚类")
        
            # 最终检查
            if current_min_dim < required_min_dim:
                raise ValueError(f"无法满足LSI聚类要求（最小维度需>{required_min_dim-1}，当前{current_min_dim}）")
    
        # 调整LSI维度范围（必须与函数默认的lsi_dim兼容）
        # 函数默认lsi_dim是list(range(1,50))，即1-49，共49个维度
        lsi_dim = list(range(1, min(49, current_min_dim-1) + 1))  # 确保不超过实际可用维度
        print(f"调整后LSI维度范围：{lsi_dim}")
    
        # 提取其他LSI参数
        lsi_params = {
            'binary': getattr(redeemr, 'lsi_binary', True),
            'res': getattr(redeemr, 'lsi_res', 0.6),
            'lsi_dim': lsi_dim,
            'rm_variants': getattr(redeemr, 'lsi_rm_variants', ["Variants310TC", "Variants3109TC", "Variants5764CT"])
            }
    
        # 运行LSI聚类（使用调整后的数据和参数）
        rep_redeemr = rep_redeemr.seurat_lsi_clustering(** lsi_params)
    except Exception as e:
        raise RuntimeError(f"代表细胞LSI聚类失败：{str(e)}")
    
    # 基于筛选后的细胞重新计算距离矩阵（确保维度匹配）
    print("\n基于筛选后的代表细胞重新计算距离矩阵...")
    try:
        cell_pct = generate_cell_pct(rep_redeemr)
        # 重新计算距离矩阵（此时基于79个代表细胞，维度应为79×79）
        rep_redeemr = add_dist(rep_redeemr, weight_df=cell_pct, lsi_dist=True)
        
        # 验证距离矩阵维度
        w_jaccard = getattr(rep_redeemr.dist_objects, 'w_jaccard', None)
        if w_jaccard is None:
            raise ValueError("add_dist未生成w_jaccard距离矩阵")
        if w_jaccard.shape[0] != len(valid_reps) or w_jaccard.shape[1] != len(valid_reps):
            raise ValueError(f"距离矩阵维度错误：{w_jaccard.shape}，应为({len(valid_reps)},{len(valid_reps)})")
        print(f"距离矩阵验证通过：{w_jaccard.shape}（与代表细胞数量匹配）")
    except Exception as e:
        raise RuntimeError(f"计算代表细胞距离矩阵失败：{str(e)}")
    
    # 构建树（此时距离矩阵维度与细胞数匹配）
    try:
        rep_redeemr = rep_redeemr.make_tree(d="w_jaccard", algorithm="nj")
    except Exception as e:
        raise RuntimeError(f"构建树失败：{str(e)}")
    
    # 可视化树（用簇ID标注）
    tree = rep_redeemr.tree.phylo
    for node in tree.find_clades():
        if node.name in valid_reps:
            # 从有效代表细胞反向找簇ID
            cluster_id = [k for k, v in representatives.items() if normalize_id(v) == normalize_id(node.name)][0]
            node.name = f"簇{cluster_id}"
    
    plt.rcParams["font.size"] = 6  # 全局字体大小（节点名称）
    plt.rcParams["axes.labelsize"] = 6  # 轴标签大小
    plt.rcParams["xtick.labelsize"] = 6
    plt.rcParams["ytick.labelsize"] = 6
    plt.figure(figsize=(150, 150))
    def format_branch_length(clade):
        """自定义分支长度格式，减小字体并简化显示"""
        if clade.branch_length:
            return f"{clade.branch_length:.3f}"  # 保留3位小数，缩短长度
        return ""
    draw(
        tree,
        branch_labels=lambda clade: f"{clade.branch_length:.4f}" if clade.branch_length else "",
        show_confidence=False,
        do_show=False  # 不自动显示，后续手动保存
    )
    plt.title("基于簇代表细胞的系统发育树", fontsize=12, pad=50)
    plt.savefig(os.path.join(output_dir, "cluster_representative_tree.png"), dpi=600, bbox_inches='tight')
    plt.close()
    print(f"系统发育树已保存至：{os.path.join(output_dir, 'cluster_representative_tree.png')}")
    
    # 绘制圆形发育树，添加布局参数
    fig, ax = draw_circular_tree(
        tree,
        label_func=lambda c: c.name if c.name else "",  # 显示节点名称
        branch_label_func=lambda clade: f"{clade.branch_length:.4f}" if clade.branch_length else ""
    )
    ax.set_title("基于簇代表细胞的圆形系统发育树", fontsize=12, pad=50)
    fig.savefig(os.path.join(output_dir, "cluster_representative_circular_tree.png"), dpi=500, bbox_inches='tight')
    plt.close()
    print(f"圆形系统发育树已保存至：{os.path.join(output_dir, 'cluster_representative_circular_tree.png')}")
    
    return rep_redeemr

# 执行分析
output_dir = "/data/project/r2py/young2_mitoConsensus/src/result/cluster_representative_analysis"
os.makedirs(output_dir, exist_ok=True)

# 传入representatives参数（用于树标注）
rep_tree_redeemr = build_tree_with_representatives(
    redeemr, 
    representative_cells, 
    representatives,  # 新增：传入簇-代表映射
    output_dir
)

# 保存结果
with open(os.path.join(output_dir, "representative_tree_result.pkl"), 'wb') as f:
    pickle.dump(rep_tree_redeemr, f)
print(f"分析结果已保存至：{output_dir}")