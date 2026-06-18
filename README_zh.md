# SeekArcTools
由 [SeekGene（寻因生物）](https://www.seekgene.com/) 开发的 SeekArcTools 是一款专为处理 SeekOne DD 单细胞多组学（ATAC + RNA）试剂盒数据而定制的生物信息学软件套装。它独立处理转录组和 ATAC 文库，用于执行诸如 Barcode 提取与校正、序列比对、Barcode 定量、联合细胞鉴定，以及生成细胞表达矩阵等任务。这使得后续能够进行聚类、差异分析，以及计算基因与 Peak 之间的关联以供下游解读。

## 安装
```
git clone https://github.com/seekgene/SeekArcTools.git
cd SeekArcTools
conda env create -n seekarctools -f env.yaml
conda activate seekarctools
pip install .
```

## 使用教程
示例 1：基本用法

设置分析所需的配置文件，包括样本数据路径、参考基因组路径等。使用以下命令运行 SeekArcTools：
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

示例 2：调整阈值以重新运行 Peak 鉴定或细胞鉴定

如果默认参数下鉴定出的 Peak 或细胞不符合要求，用户可以调整参数跳过前面的步骤（如比对）并重新运行流程，从而在优化结果的同时节省时间。
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

* 注意：确保之前运行产生的所有文件保持完整，且未被删除或移动。`--outdir` 参数需指定为最初运行 seekarctools 分析时的目录路径。`--min_atac_count` 必须与 `--min_gex_count` 配合使用；单独使用将不会生效。

## 参数说明

| 参数名            | 参数说明                                                     |
| ----------------- | ------------------------------------------------------------ |
| --rnafq1          | RNA 文库 R1 端 fastq 文件路径                                  |
| --rnafq2          | RNA 文库 R2 端 fastq 文件路径                                  |
| --atacfq1         | ATAC 文库 R1 端 fastq 文件路径                                  |
| --atacfq2         | ATAC 文库 R2 端 fastq 文件路径                                  |
| --samplename      | 样本名称                                                   |
| --chemistry       | 试剂化学类型。可用选项: DD_AG, DD5_AG                                                 |
| --outdir          | 输出目录路径。默认值: ./                                   |
| --skip_misB       | 若启用，则 Barcode 不允许存在碱基错配。默认允许 1 个错配。                |
| --skip_misL       | 若启用，则 Linker 不允许存在碱基错配。默认允许 1 个错配。                 |
| --skip_multi      | 若启用，则丢弃可被校正到白名单中多个 Barcode 的 reads。默认会校正到频率最高的 Barcode。 |
| --skip_len        | 在接头过滤后跳过对短序列的过滤，短序列将被保留使用。                 |
| --core            | 用于分析的线程数。                                           |
| --include-introns | 禁用时，仅使用外显子 reads 进行定量。启用时，内含子 reads 也会被用于定量。 |
| --refpath         | 参考基因组目录路径。                                             |
| --star_path       | 外部 STAR 软件路径。如果参考基因组中的索引是由其他版本的 STAR 构建的，请指定其路径。 |
| --qvalue          | Peak 鉴定的最小 FDR (q-value) 阈值。默认值: 0.001。               |
| --nolambda        | 若为 True，MACS3 将使用背景 lambda 作为局部 lambda。这意味着 MACS3 不会考虑 Peak 候选区域的局部偏倚。 |
| --snapshift       | MACS3 Peak 鉴定偏移量 (shift size)。默认值: 0。                         |
| --extsize         | MACS3 Peak 鉴定延伸大小 (extension size)。默认值: 400。                       |
| --min_len         | Peak 的最小长度。如果未设置，将默认等同于 “extsize”。默认值: 400。       |
| --broad           | 若启用，则执行 broad peak 鉴定并生成 UCSC gappedPeak 格式的结果，该格式包含嵌套的 Peak 结构。 |
| --broad_cutoff    | Broad peak 鉴定的阈值。默认值: 0.1。                          |
| --min_atac_count  | 细胞鉴定覆盖参数：定义一个细胞 Barcode 在 Peak 区域内最小的 ATAC 转座事件数 (ATAC counts)。       |
| --min_gex_count   | 细胞鉴定覆盖参数：定义一个细胞 Barcode 最小的基因表达 UMI 数 (GEX counts)。                   |
| -h, --help        | 显示此参数说明帮助信息。                                                 |

## 输出说明
以下为输出目录结构：每一行代表一个文件或文件夹，由 “├──” 指示，后面的数字标记了重要的输出文件。
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
7. ATAC 文库的比对 bam 文件
8. ATAC 鉴定的 Peaks 文件
9. RNA 文库的比对 bam 文件
10. rds 格式的矩阵及分析对象
11. 最终 HTML 分析报告
12. CSV 格式的质量控制统计信息

## 许可证
本项目采用 MIT 许可证 - 有关详细信息，请参阅 LICENSE 文件。