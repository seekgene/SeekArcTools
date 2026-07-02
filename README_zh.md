# SeekArcTools
由 [SeekGene（寻因生物）](https://www.seekgene.com/) 开发的 SeekArcTools 是一款专为处理 SeekOne DD 单细胞多组学（ATAC + RNA）试剂盒数据而定制的生物信息学软件套装，同时支持该体系下 RNA 文库的单独分析。它可分别处理转录组和 ATAC 文库，用于执行 Barcode 提取与校正、序列比对、Barcode 定量、细胞鉴定，以及生成表达矩阵、Peak 矩阵、Seurat RDS 和分析报告等任务。输出结果可用于后续聚类、差异分析等下游分析；对于 ATAC + RNA 多组学数据，还可进一步支持基因与 Peak 之间的关联分析。

## 安装
```
git clone https://github.com/seekgene/SeekArcTools.git
cd SeekArcTools
conda env create -f env.yaml
conda activate seekarctools
pip install .
```

## 使用教程
示例 1：ARC 多组学分析

ARC 多组学分析同时输入 RNA 文库和 ATAC 文库 FASTQ 文件，用于完成 RNA 定量、ATAC Peak 鉴定、联合细胞鉴定、矩阵生成及报告输出。

```shell
seekarctools arc run \
    --rnafq1 /path/to/demo/demo_GE_S1_L001_R1_001.fastq.gz \
    --rnafq2 /path/to/demo/demo_GE_S1_L001_R2_001.fastq.gz \
    --atacfq1 /path/to/demo/demo_arc_S2_L002_R1_001.fastq.gz \
    --atacfq2 /path/to/demo/demo_arc_S2_L002_R2_001.fastq.gz \
    --samplename demo \
    --outdir /path/to/outdir \
    --refpath /path/to/reference/GRCh38 \
    --chemistry DD_AG \
    --include-introns \
    --core 16
```

示例 2：RNA 文库单独分析

RNA 文库单独分析用于仅包含 RNA/GEX 文库的数据。流程会完成 Barcode 处理、STAR 比对、UMI 定量、细胞鉴定、表达矩阵生成、Seurat RDS 生成和 GEX 质控统计输出。

```shell
seekarctools rna run \
    --rnafq1 /path/to/demo/demo_GE_S1_L001_R1_001.fastq.gz \
    --rnafq2 /path/to/demo/demo_GE_S1_L001_R2_001.fastq.gz \
    --samplename demo \
    --outdir /path/to/outdir \
    --refpath /path/to/reference/GRCh38 \
    --chemistry DD_AG_RNA \
    --include-introns \
    --core 16
```

示例 3：调整阈值以重新运行 Peak 鉴定或细胞鉴定

如果默认参数下鉴定出的 Peak 或细胞不符合要求，可以使用 `arc retry` 跳过前面的步骤（如比对），基于已有中间结果重新运行后续流程，从而在优化结果的同时节省时间。

```shell
seekarctools arc retry \
    --samplename demo \
    --outdir /path/to/outdir \
    --refpath /path/to/reference/GRCh38 \
    --core 16 \
    --qvalue 0.01 \
    --snapshift 0 \
    --extsize 200 \
    --min_len 200 \
    --min_atac_count 1000 \
    --min_gex_count 500
```

* 注意：确保之前运行产生的所有文件保持完整，且未被删除或移动。`--outdir` 参数需指定为最初运行  `seekarctools arc run` 时的输出目录路径。`--min_atac_count` 必须与 `--min_gex_count` 配合使用；单独使用将不会生效。

## 参数说明

| 参数名                     | 参数说明                                                        |
| ----------------------- | ----------------------------------------------------------- |
| `--rnafq1`              | RNA 文库 R1 端 FASTQ 文件路径                                      |
| `--rnafq2`              | RNA 文库 R2 端 FASTQ 文件路径                                      |
| `--atacfq1`             | ATAC 文库 R1 端 FASTQ 文件路径                                     |
| `--atacfq2`             | ATAC 文库 R2 端 FASTQ 文件路径                                     |
| `--samplename`          | 样本名称                                                        |
| `--chemistry`           | 试剂化学类型。ARC 多组学分析使用 `DD_AG` 或 `DD5_AG`；RNA 文库单独分析使用 `DD_AG_RNA` |
| `--outdir`              | 输出目录路径。默认值：`./`                                             |
| `--skip_misB`           | 启用后，Barcode 不允许存在碱基错配；默认允许 1 个错配                            |
| `--skip_misL`           | 启用后，Linker 不允许存在碱基错配；默认允许 1 个错配                             |
| `--skip_multi`          | 启用后，丢弃可被校正到多个白名单 Barcode 的 reads；默认校正到频率最高的 Barcode         |
| `--skip_len`            | 启用后，跳过接头过滤后的短序列过滤，保留短序列用于后续分析                               |
| `--core`                | 用于分析的线程数                                                    |
| `--memory`              | RNA 文库单独分析中生成 Seurat RDS 使用的内存参数                            |
| `--include-introns`     | 启用后，内含子 reads 也会被用于基因表达定量；默认仅使用外显子 reads                    |
| `--refpath`             | 参考基因组目录路径                                                   |
| `--star_path`           | 外部 STAR 软件路径。如果参考基因组索引由指定 STAR 版本构建，可指定该路径                  |
| `--scoremin`            | STAR 宽松比对参数，对应 `--outFilterScoreMinOverLread`               |
| `--matchnmin`           | STAR 宽松比对参数，对应 `--outFilterMatchNminOverLread`              |
| `--pairendAlignment`    | RNA 文库单独分析中启用 STAR paired-end 比对；默认使用 R2 单端比对               |
| `--pairendMinR1Len`     | RNA paired-end 比对时保留 R1 的最小长度。默认值：`20`                      |
| `--expect_cells`        | RNA 文库单独分析中预期回收细胞数。默认值：`3000`                               |
| `--pvalue`              | RNA 细胞鉴定中 ambient RNA 检验的 p-value 阈值。默认值：`0.01`             |
| `--force_cells`         | RNA 文库单独分析中强制取 UMI 数最高的前 N 个 Barcode 作为细胞；默认值：`0`，表示自动识别    |
| `--qvalue`              | Peak 鉴定的最小 FDR（q-value）阈值。默认值：`0.001`                       |
| `--nolambda`            | 启用后，MACS3 使用背景 lambda 作为 local lambda，不考虑 Peak 候选区域的局部偏倚    |
| `--snapshift`           | MACS3 Peak 鉴定偏移量（shift size）。默认值：`0`                        |
| `--extsize`             | MACS3 Peak 鉴定延伸大小（extension size）。默认值：`400`                 |
| `--min_len`             | Peak 的最小长度。如果未设置，将默认等同于 `extsize`                           |
| `--broad`               | 启用后，执行 broad peak 鉴定并生成 UCSC gappedPeak 格式结果                |
| `--broad_cutoff`        | Broad peak 鉴定阈值。默认值：`0.1`                                   |
| `--min_atac_count`      | 细胞鉴定覆盖参数：定义细胞 Barcode 在 Peak 区域内最小的 ATAC 转座事件数（ATAC counts） |
| `--min_gex_count`       | 细胞鉴定覆盖参数：定义细胞 Barcode 最小的基因表达 UMI 数（GEX counts）             |
| `--keep_mito`           | 启用后，保留线粒体 reads；默认过滤线粒体 reads                               |
| `--strict_cell_calling` | 启用后，在 K-Means 细胞鉴定基础上增加 ambient RNA 背景检验，以获得更严格的细胞识别结果      |
| `--max_cells`           | 限制自动识别细胞数上限；未设置时不限制                                         |
| `-h`, `--help`          | 显示命令帮助信息                                                    |

## 输出说明

### ARC 多组学分析输出

以下为 ARC 多组学分析的主要输出目录结构。每一行代表一个文件或文件夹，由 “├──” 指示，后面的数字标记了重要输出文件。

```shell
./
├── outs
│   ├── filtered_feature_bc_matrix              1
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── filtered_peaks_bc_matrix                2
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── raw_feature_bc_matrix                   3
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── raw_peaks_bc_matrix                     4
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── demo_A_fragments.tsv.gz                 5
│   ├── demo_A_fragments.tsv.gz.tbi             6
│   ├── demo_A_mem_pe_Sort.bam                  7
│   ├── demo_A_peaks.bed                        8
│   ├── demo_E_SortedByName.bam                 9
│   ├── demo.rds                                10
│   ├── demo_report.html                        11
│   └── demo_summary.csv                        12
└── analysis

```
1. 转录组细胞过滤后稀疏矩阵
2. ATAC 细胞过滤后稀疏矩阵
3. 转录组原始稀疏矩阵
4. ATAC 原始稀疏矩阵
5. ATAC 每个片段(fragment)信息文件
6. Fragments 文件索引
7. ATAC 文库比对 BAM 文件
8. ATAC Peaks 文件
9. RNA 文库比对 BAM 文件
10. RDS 格式的 Seurat 对象及下游分析结果
11. 最终 HTML 分析报告
12. CSV 格式的质量控制统计信息

### RNA 文库单独分析输出

以下为 RNA 文库单独分析的主要输出目录结构。

```shell
./
├── analysis
│   └── demo_E
│       ├── step1
│       │   ├── demo_E_1.fq.gz
│       │   └── demo_E_2.fq.gz
│       ├── step2
│       │   └── featureCounts
│       │       └── demo_E_SortedByName.bam     1
│       ├── step3
│       │   ├── raw_feature_bc_matrix           2
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   ├── filtered_feature_bc_matrix      3
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   ├── counts.xls
│       │   ├── detail.xls
│       │   └── umi.xls
│       └── step4
│           ├── demo.rds                        4
│           ├── gex_tsne_umi.xls                5
│           └── gex_FindAllMarkers.xls          6
└── outs
    └── demo_summary.csv                        7
```

1. RNA 文库比对 BAM 文件
2. RNA 原始表达矩阵
3. RNA 过滤后细胞表达矩阵
4. RDS 格式的 Seurat 对象
5. 细胞降维坐标及 UMI/线粒体比例等信息
6. 基于 Seurat 聚类结果的 marker 基因表
7. RNA/GEX 质控统计信息


## 许可证
本项目采用 MIT 许可证 - 有关详细信息，请参阅 LICENSE 文件。
