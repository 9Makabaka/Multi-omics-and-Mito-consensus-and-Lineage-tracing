from io import StringIO
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from Bio import Phylo

def draw_circular_tree(tree, label_func=str, branch_label_func=None, figsize=(15, 15)):
    """
    绘制圆形布局的系统发育树
    
    参数:
        tree: 系统发育树对象（如Bio.Phylo的Tree对象）
        label_func: 节点标签处理函数
        branch_label_func: 分支标签处理函数（输入clade，返回标签文本，可选）
        figsize: 图像大小
    """
    # 计算节点坐标（极坐标：半径和角度）
    def get_radii_and_angles(tree):
        depths = tree.depths()
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        
        terminals = tree.get_terminals()
        n_terminals = len(terminals)
        angles = {}
        for i, tip in enumerate(terminals):
            angles[tip] = 2 * np.pi * (i + 1) / (n_terminals + 1)
        
        def calc_internal_angle(clade):
            for child in clade:
                if child not in angles:
                    calc_internal_angle(child)
            angles[clade] = np.mean([angles[child] for child in clade.clades])
        
        if tree.root.clades:
            calc_internal_angle(tree.root)
        
        return depths, angles
    
    radii, angles = get_radii_and_angles(tree)
    
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_xticks([])
    ax.set_yticks([])
    
    lines = []  
    
    def draw_clade(clade, parent_r, parent_theta):
        r = radii[clade]
        theta = angles[clade]
        
        # 绘制径向线
        lines.append([(parent_theta, parent_r), (theta, r)])
        
        # 绘制子节点间的弧线
        if clade.clades:
            child_thetas = [angles[child] for child in clade.clades]
            for i in range(len(child_thetas)-1):
                start_theta = child_thetas[i]
                end_theta = child_thetas[i+1]
                arc_thetas = np.linspace(start_theta, end_theta, 20)
                arc_rs = [r] * 20
                arc_points = list(zip(arc_thetas, arc_rs))
                lines.append(arc_points)
        
        # 节点标签
        label = label_func(clade)
        if label:
            label_r = r * 1.05
            ax.text(theta, label_r, label, ha='center', va='center', fontsize=8)
        
        # 分支标签（使用传入的处理函数）
        if branch_label_func is not None and clade.branch_length is not None:
            mid_theta = (parent_theta + theta) / 2
            mid_r = (parent_r + r) / 2
            # 调用分支标签函数处理
            ax.text(mid_theta, mid_r, branch_label_func(clade),
                    ha='center', va='center', fontsize=6)
        
        # 递归绘制子节点
        for child in clade:
            draw_clade(child, r, theta)
    
    draw_clade(tree.root, 0, 0)
    ax.add_collection(LineCollection(lines, color='black', linewidths=1))
    
    max_r = max(radii.values())
    ax.set_ylim(0, max_r * 1.1)
    
    return fig, ax

# 示例使用
if __name__ == "__main__":
    output_dir = "/data/project/r2py/young2_mitoConsensus/src/result"
    os.makedirs(output_dir, exist_ok=True)
    # 加载示例树（Newick格式）
    newick = "(A:0.1, (B:0.2, C:0.3):0.4, (D:0.5, (E:0.6, F:0.7):0.8):0.9);"
    tree = Phylo.read(StringIO(newick), "newick")
    
    # 绘制圆形树
    fig, ax = draw_circular_tree(
        tree,
        label_func=lambda c: c.name if c.name else "",  # 显示节点名称
        branch_label_func=lambda clade: f"{clade.branch_length:.4f}" if clade.branch_length else ""
    )
    plt.title("圆形系统发育树示例", y=1.1)  # 标题
    plt.savefig(os.path.join(output_dir, "example_tree.png"), dpi=600, bbox_inches='tight')
    plt.close()
    print(f"示例发育树已保存至：{os.path.join(output_dir, 'example_tree..png')}")
    