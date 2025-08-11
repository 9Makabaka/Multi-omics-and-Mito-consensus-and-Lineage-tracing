import os
import pickle
from matplotlib import pyplot as plt
from mito_consensus import create_redeemr,read_rds_variants,make_matrix,generate_cell_pct
from dist_cal import add_dist, from_dist2graph
import igraph as ig
import scanpy as sc
from Bio.Phylo import draw

from visual import draw_circular_tree

if __name__ == "__main__":
    # # 输入文件路径
    # rds_file_path = "/data/dataset/young2_HSC/Youn2.HSC.Consensus.final/VariantsGTSummary.RDS"
    
    # # 步骤1：读取RDS文件（提取"Sensitive"严格度的数据）
    # print(f"读取文件: {rds_file_path}，提取严格度：S（Sensitive）")
    # variants_gtsummary = read_rds_variants(
    #     rds_path=rds_file_path,
    #     thr="S"
    # )
    
    # # 步骤2：创建redeemR对象
    # redeemr = create_redeemr(
    #     variants_gtsummary=variants_gtsummary,
    #     qualified_cell_cut=10,
    #     vaf_cut=1,
    #     cell_cut=2,
    #     maxcts_cut=2,
    #     thr="S"
    # )
    
    # # 步骤3：生成细胞-变异矩阵
    # redeemr = make_matrix(redeemr, onlyhetero=True)
    
    # # 步骤4：执行LSI聚类分析
    # print("\n执行LSI聚类...")
    # redeemr = redeemr.seurat_lsi_clustering(
    #     binary=True,
    #     res=0.6,
    #     lsi_dim=list(range(1, 50)),  # 使用第2到第50个LSI维度（对应R的2:50）
    #     rm_variants=["Variants310TC", "Variants3109TC", "Variants5764CT"]
    # )
    # adata = redeemr.seurat  # 获取AnnData对象
    # # 绘制UMAP聚类图
    # plt.figure(figsize=(8, 6))
    # sc.pl.umap(
    #     adata,
    #     color="leiden",  # 按聚类标签着色（leiden聚类结果）
    #     title="LSI Clustering (UMAP)",
    #     show=False  # 关闭自动显示，后续统一控制
    # )
    # plt.savefig("/data/project/r2py/young2_mitoConsensus/src/result/lsi_clustering_umap.png", dpi=300, bbox_inches="tight")  # 保存图片
    # plt.show()  # 显示图片
    
    # # 步骤5：生成CellPCT（权重数据）
    # print("\n生成CellPCT权重数据...")
    # cell_pct = generate_cell_pct(redeemr)
    
    # print("检查LSI结果是否存在：")
    # print("redeemr.seurat是否为None？", redeemr.seurat is None)
    # if redeemr.seurat is not None:
    #     print("LSI降维数据是否存在？", 'X_pca_subset' in redeemr.seurat.obsm)
    #     print("Leiden聚类结果是否存在？", 'leiden' in redeemr.seurat.obs)
    
    # # 步骤6：计算并添加距离指标
    # print("\n计算距离指标...")
    # redeemr = add_dist(
    #     redeemr_obj=redeemr,
    #     weight_df=cell_pct,  # 使用自动生成的权重数据
    #     lsi_dist=True  # 基于之前的LSI结果计算距离
    # )
    # 保存redeemr对象
    save_path_dist = "/data/project/r2py/young2_mitoConsensus/src/result/redeemr_with_dist.pkl"
    # os.makedirs(os.path.dirname(save_path_dist), exist_ok=True)
    # with open(save_path_dist, 'wb') as f:
    #     pickle.dump(redeemr, f)
    # print(f"redeemr对象已保存至：{save_path_dist}")
    
    # 加载保存的redeemr对象
    with open(save_path_dist, 'rb') as f:
        redeemr = pickle.load(f)
    print("redeemr对象加载完成")
    
    # 步骤7：生成互近邻（MNN）网络
    print("\n生成互近邻（MNN）网络...")
    distance_matrix = redeemr.dist_objects.w_jaccard
    mnn_graph = from_dist2graph(
        d=distance_matrix,
        k_param=30,  # 每个细胞找30个最近邻
        return_igraph=True
    )
    # # 可视化网络
    # layout = mnn_graph.layout_kamada_kawai()  # 计算布局（网络节点的空间排列）
    # # 绘图
    # ig.plot(
    #     mnn_graph, 
    #     layout=layout, 
    #     bbox=(800, 800),  # 图像大小
    #     vertex_size=10,   # 节点大小
    #     vertex_label_size=8,  # 节点标签大小
    #     target="/data/project/r2py/young2_mitoConsensus/src/result/mnn_network.png"  # 保存为图片文件
    # )
    
    # 步骤8：系统发育树构建
    print("\n系统发育树构建...")
    redeemr = redeemr.make_tree(d="w_jaccard", algorithm="nj")
    # 关键步骤：保存redeemr对象（类似R的saveRDS）
    save_path = "/data/project/r2py/young2_mitoConsensus/src/result/redeemr_with_tree.pkl"
    # 确保保存目录存在（不存在则创建）
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    # 序列化并保存对象
    with open(save_path, 'wb') as f:
        pickle.dump(redeemr, f)
    print(f"redeemr对象已保存至：{save_path}")
    
    # （后续使用时）加载保存的redeemr对象（类似R的readRDS）
    # with open(save_path, 'rb') as f:
    #     redeemr = pickle.load(f)
    # print("redeemr对象加载完成")
    
    tree = redeemr.tree.phylo
    # 绘制系统发育树
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
    plt.title("系统发育树", fontsize=12, pad=50)
    plt.savefig("/data/project/r2py/young2_mitoConsensus/src/result/phylogenetic_tree.png", dpi=600, bbox_inches='tight')
    plt.close()
    print(f"系统发育树已保存至：/data/project/r2py/young2_mitoConsensus/src/result/phylogenetic_tree.png")
    
    # 绘制圆形发育树，添加布局参数
    fig, ax = draw_circular_tree(
        tree,
        label_func=lambda c: c.name if c.name else "",  # 显示节点名称
        branch_label_func=lambda clade: f"{clade.branch_length:.4f}" if clade.branch_length else ""
    )
    ax.set_title("圆形系统发育树", fontsize=12, pad=50)
    fig.savefig("/data/project/r2py/young2_mitoConsensus/src/result/cluster_representative_circular_tree.png", dpi=500, bbox_inches='tight')
    plt.close()
    print(f"圆形系统发育树已保存至：/data/project/r2py/young2_mitoConsensus/src/result/cluster_representative_circular_tree.png")

        
    # 查看结果
    print("\n===== 基础数据概览 =====")
    print("合格细胞元数据:")
    print(redeemr.cell_meta.head())
    print("\n过滤后的变异特征:")
    print(redeemr.v_filtered.head())
    print("\n纯质性变异（前5个）:")
    print(redeemr.homo_variants[:5] if redeemr.homo_variants else "无")
    
    print("\n===== 矩阵信息 =====")
    if redeemr.cts_mtx:
        print(f"计数矩阵形状: {redeemr.cts_mtx['matrix'].shape}")
        print(f"前3个细胞: {redeemr.cts_mtx['rows'][:3]}")
        print(f"前3个变异: {redeemr.cts_mtx['columns'][:3]}")
    
    if redeemr.cts_mtx_bi:
        print(f"二值化矩阵形状: {redeemr.cts_mtx_bi['matrix'].shape}")
        
        # 查看聚类结果
    print("\n===== 聚类结果概览 =====")
    if redeemr.seurat is not None:
        print("聚类数量分布:")
        print(redeemr.seurat.obs["leiden_clusters"].value_counts())
        print("\nUMAP坐标前5行:")
        print(redeemr.seurat.obsm["X_umap"][:5, :])
    
    # 查看距离结果
    print("\n===== 距离计算结果 =====")
    if redeemr.dist_objects:
        print(f"Jaccard距离矩阵形状: {redeemr.dist_objects.jaccard.shape if redeemr.dist_objects.jaccard is not None else '未计算'}")
        print(f"加权Jaccard距离矩阵形状: {redeemr.dist_objects.w_jaccard.shape if redeemr.dist_objects.w_jaccard is not None else '未计算'}")
        print(f"LSI距离矩阵形状: {redeemr.dist_objects.lsi_dist.shape if redeemr.dist_objects.lsi_dist is not None else '未计算'}")
        
    # 查看网络基本信息
    print("\n===== MNN网络信息 =====")
    if mnn_graph:
        print(f"细胞数量: {mnn_graph.vcount()}")
        print(f"互近邻边数量: {mnn_graph.ecount()}")

