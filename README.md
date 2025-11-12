# SeekArcTools
SeekArcTools is a software developed by SEEKGENE for processing single-cell ATAC multi-omics data. The software independently processes transcriptome libraries and ATAC libraries for tasks, such as:Barcode extraction and correction, alignment.Quantification of each barcode, joint cell identification, and generation of cell expression matrices for downstream analysis.Enables subsequent clustering and differential analysis, as well as calculation of associations between genes and peaks for downstream interpretation.
## Installation
```
git clone https://github.com/seekgene/SeekArcTools.git
cd SeekArcTools
conda env create -n seekarctools -f env.yaml
conda activate seekarctools
pip install .
```

## Tutorials
Example 1: Basic usage

Set up the necessary configuration files for the analysis, including sample data paths, reference genome paths, etc. Run the SeekArcTools using the following command:
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

Example 2: Adjusting thresholds for re-running peak calling or cell calling

If the identified peaks or cells under default parameters do not meet requirements, users can adjust parameters to skip preceding steps (e.g., alignment) and re-run the workflow to optimize results while saving time.
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

* Note: Ensure all files from previous runs remain intact and are not deleted or relocated. The `--outdir` specifies the directory path from the initial seekarctools analysis. The `--min_atac_count` must be used together with `--min_gex_count`; it will not take effect if used alone.

## Parameter descriptions

| Parameters        | Parameter descriptions                                                     |
| ----------------- | ------------------------------------------------------------ |
| --rnafq1          | Paths to R1 fastq files of RNA library                                  |
| --rnafq2          | Paths to R2 fastq files of RNA library                                  |
| --atacfq1         | Paths to R1 fastq files of ATAC library                                  |
| --atacfq2         | Paths to R2 fastq files of ATAC library                                  |
| --samplename      | Sample name                                                   |
| --chemistry       | Reagent type. Available options: DD_AG, DD5_AG                                                 |
| --outdir          | output directory. Default: ./                                   |
| --skip_misB       | If enabled, no base mismatch is allowed for barcode. Default is 1.                |
| --skip_misL       | If enabled, no base mismatch is allowed for linker. Default is 1.                 |
| --skip_multi      | If enabled, discard reads that can be corrected to multiple white-listed barcodes. Barcodes are corrected to the barcode with the highest frequency by default. |
| --skip_len        | Skip filtering short reads after adapter filter, short reads will be used.                 |
| --core            | Number of threads used for the analysis.                                           |
| --include-introns | When disabled, only exon reads are used for quantification. When enabled, intron reads are also used for quantification. |
| --refpath         | The path of reference genome.                                             |
| --star_path       | External STAR software path. If the index in the reference genome is built by another STAR, please specify its path. |
| --qvalue          | Minimum FDR (q-value) cutoff for peak detection. Default: 0.001.               |
| --nolambda        | If True, MACS3 will use the background lambda as local lambda. This means MACS3 will not consider the local bias at peak candidate regions. |
| --snapshift       | MACS3 peak detection shift size. Default: 0.                         |
| --extsize         | MACS3 peak detection extension size. Default: 400.                       |
| --min_len         | Minimum length of peaks. If not set, it will be set to “extsize”. Default: 400.       |
| --broad           | If enabled, perform broad peak calling and generate results in UCSC gappedPeak format, which encapsulates the nested structure of peaks. |
| --broad_cutoff    | Threshold for broad peak calling. Default: 0.1.                          |
| --min_atac_count  | Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.       |
| --min_gex_count   | Cell caller override: define the minimum number of GEX UMI counts for a cell barcode.                   |
| -h, --help        | Show this parameter descriptions.                                                 |

## Output descriptions
“Here’s the output directory structure: each line represents a file or folder, indicated by “├──”, and the numbers indicate three important output files.
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
1. Transcriptome cell sparse matrix
2. ATAC cell sparse matrix
3. Transcriptome raw sparse matrix
4. ATAC raw sparse matrix
5. ATAC per fragment information file
6. Fragments file index
7. The bam file of ATAC libiary
8. ATAC peaks
9. The bam file of RNA libiary
10. Matrix in rds format
11. Final report in html
12. Quality control information in csv

## License
This project is licensed under the MIT License - see the LICENSE file for details.