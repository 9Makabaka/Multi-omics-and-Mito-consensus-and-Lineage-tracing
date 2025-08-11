线粒体DNA（mtDNA）的共识变异检测（consensus variant calling）：

1. **单细胞级分析**：基于细胞barcode（如10X Genomics的`BC`标签）区分不同细胞。
2. **分子级统计**：将读段按分子（`CellBC_Start_End`）分组，避免PCR重复干扰。
3. **多级过滤**：输出4种严格级别的变异结果（Total/VerySensitive/Sensitive/Specific）。
4. **链特异性支持**：记录正反链覆盖信息，提高变异检测可靠性。

### 1.1输入输出

**输入文件**：

- `barcodes_file`：细胞条形码列表（每个条形码对应一个细胞）；
- `bam_file`：比对后的 BAM 文件（存储测序读段的序列、质量、位置等信息）；
- `mito_ref_file`：线粒体参考基因组（用于判断变异是否与参考序列不同）。

**输出文件**：

- 4 个变异结果文件（按严格度区分：Total、VerySensitive、Sensitive、Specific）；
- 1 个计数文件（`QualifiedTotalCts`，记录每个细胞每个位置的支持数）。

### 1.2核心数据结构

**关键字典**：通过解析 BAM 文件，构建两个核心字典，用于关联测序读段与对应的 DNA 分子。

- **ReadPairDict**：键为读段对名称（read pair name），值为该读段对包含的两个读段（read1 和 read2）。作用是快速通过读段对名称获取具体的测序序列、质量、比对位置等信息。
- **MoleculeDict**：键为 “分子 ID”（由细胞条形码 + 起始位置 + 结束位置组成，格式：`CellBC_pos_start_pos_end`），值为该分子对应的所有读段对名称。

```python
ReadPairDict = {read_name: [read1, read2]}  # 读段名称到双端读段的映射MoleculeDict = {CellBC_Start_End: [read_pair1, ...]}  # 分子到读段对的映射
```

> 这里的 “分子” 指的是同一个 DNA 分子经测序得到的所有读段集合（同一分子的读段应覆盖相同区域）。构建逻辑：
> 
> 1. 过滤无效读段对（非成对、链方向相同的读段对）；
> 2. 基于读段的模板长度（tlen）计算分子的起始和结束位置；
> 3. 将同一分子的读段对关联到对应的分子 ID。

**统计矩阵**

| 矩阵名称 | 维度 | 说明 |
| --- | --- | --- |
| `SG_Genotypes` | `(max_bp, 5)` | 单链支持的碱基数（A/C/G/T/N） |
| `DB_Genotypes` | `(max_bp, 5)` | 双链支持的碱基数 |
| `Strand_mtx` | `(max_bp, 2)` | 正/反链支持数（0:+，1:-） |
| `TotalMoleculeCtsMatrix` | `(n_cells, max_bp, 4)` | 各细胞各位置的覆盖深度（4种过滤级别） |

### 1.3核心算法

**1.3.1数据预处理**

对分子包含的每个读段对（`read_pair`），解析其序列、质量、比对位置，区分**重叠区域**和**非重叠区域**，分别更新单链和双链矩阵：

- **读段对解析**：
    - `seq_0`/`seq_1`：读段 1 和读段 2 的序列；
    - `quality_0`/`quality_1`：读段 1 和读段 2 的碱基质量值；
    - `pos_array_0`/`pos_array_1`：读段 1 和读段 2 与参考基因组比对的位置数组（仅匹配的位置）。
- **区域划分**：
    - `pos_array_overlap`：读段 1 和读段 2 重叠的碱基位置；
    - `pos_array_specific_0`/`pos_array_specific_1`：读段 1 和读段 2 的非重叠区域。

**1.3.2碱基统计**

```python
for read_pair in MoleculeDict[m]:
    # 分离单链和双链覆盖区域    pos_array_overlap = 重叠区域
    pos_array_specific_0 = 读段0特有区域    pos_array_specific_1 = 读段1特有区域    # 统计单链支持（仅高质量碱基）    if quality > BaseQ_thld_hi:
        SG_Genotypes[pos, base] += 1    # 统计双链支持（需双链一致）    if seq0 == seq1 and (quality0 > threshold or quality1 > threshold):
        DB_Genotypes[pos, base] += 1
```

- **填充单链矩阵（SG_Genotypes）**：
    - 非重叠区域的碱基：若碱基质量 > 阈值（`BaseQ_thld_hi`），则在对应位置的碱基列 + 1；否则记为 N（第 4 列）；
    - 同时更新链方向矩阵（`Strand_mtx`），记录该碱基来自正向还是反向链。
- **填充双链矩阵（DB_Genotypes）**：
    - 重叠区域的碱基：仅当两条读段的碱基一致，且至少一条读段的碱基质量 > 阈值时，在对应位置的碱基列 + 1；否则记为 N；
    - 若碱基不一致，直接记为 N（代码中注释了冲突处理逻辑，可能因可靠性低而舍弃）。

**1.3.3变异检测**

对每个碱基位置（`i`），基于`SG_Genotypes`和`DB_Genotypes`的总和，计算基因型并判断是否为变异。

```python
Cur_Genotype_array = SG + DB  # 合并单双链支持FamSize = sum(Cur_Genotype_array)  # 总支持数Call = argmax(Cur_Genotype_array)   # 最可能碱基CSS = GT_Cts / FamSize             # 共识支持分数# 多级过滤逻辑if DB_Cts > 0:  # 双链支持    VerySensitive: CSS > 0.75 & FamSize >= 1    Specific: CSS > 0.9 & FamSize >= 3else:           # 仅单链支持    VerySensitive: CSS > 0.75 & FamSize >= 2    Specific: CSS > 0.9 & FamSize >= 4
```

**1.3.4输出结果**

- **变异记录格式**：
    
    ```
    Molecule CellBC Position Variant Call Ref FamSize GT_Cts CSS DB_Cts SG_Cts StrandF StrandR
    ```
    
- **覆盖深度文件**：
    
    ```
    CellBC Position TotalDepth VerySensitiveDepth SensitiveDepth SpecificDepth
    ```