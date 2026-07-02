suppressPackageStartupMessages({
  library(future)
  library(Signac)
  library(Seurat)
  library(dplyr)
  library(argparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(Rsamtools)
})

parser = ArgumentParser()
parser$add_argument("--gex_matrix", help="gex. step3/filtered_feature_bc_matrix")
parser$add_argument("--atac_matrix", help="atac. step3/filtered_peaks_bc_matrix")
parser$add_argument("--fragpath", help="atac. step3/asample_fragments.tsv.gz")
parser$add_argument("--rawname", help="raw name")
parser$add_argument("--samplename", help="sample name")
parser$add_argument("--outdir", help="outdir")
parser$add_argument("--core", help="Parallel running cores")
parser$add_argument("--memory", help="Memory usage")
parser$add_argument("--ref_gtf", help="reference. gene/gtf")
parser$add_argument("--ref_fa", help="reference. fasta/genome.fa")
args <- parser$parse_args()

outdir=args$outdir
gex_matrix=args$gex_matrix
atac_matrix=args$atac_matrix
fragpath=args$fragpath
rawname=args$rawname
samplename=args$samplename
core=args$core
memory=args$memory
ref_gtf=args$ref_gtf
ref_fa=args$ref_fa

cat("available memory:", memory, "\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# multicore
options(future.globals.maxSize = as.integer(memory) * 1024 ^ 3)
options(future.rng.onMisuse = "ignore")
plan("multicore", workers = as.integer(core))
plan()

gex_data <- Read10X(data.dir = gex_matrix)
obj <- CreateSeuratObject(counts = gex_data, assay = "RNA")
atac_data <- Read10X(data.dir = atac_matrix)

# read gtf using fast and robust method (data.table + stringi)
cat("Reading GTF using fast and robust method...\n")
if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("stringi", quietly = TRUE)) {
  cat("data.table or stringi not found, falling back to robust read.table method...\n")
  gtf_df <- read.table(ref_gtf, sep = "\t", header = FALSE, quote = "", comment.char = "#", stringsAsFactors = FALSE)
  colnames(gtf_df) <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  
  extract_attr <- function(text, key) {
    pattern <- paste0(key, ' "([^"]*)"')
    matches <- regexec(pattern, text)
    res <- regmatches(text, matches)
    sapply(res, function(x) if(length(x) > 1) x[2] else NA)
  }
  
  gtf_df$gene_id <- extract_attr(gtf_df$attributes, "gene_id")
  gtf_df$gene_name <- extract_attr(gtf_df$attributes, "gene_name")
  gtf_df$transcript_id <- extract_attr(gtf_df$attributes, "transcript_id")
  gtf_df$gene_biotype <- extract_attr(gtf_df$attributes, "gene_biotype")
  if (all(is.na(gtf_df$gene_biotype))) {
    gtf_df$gene_biotype <- extract_attr(gtf_df$attributes, "gene_type")
  }
} else {
  library(data.table)
  library(stringi)
  # 使用 cmd 直接在系统层面过滤掉无用 feature，极大节省内存
  gtf_df <- fread(
    cmd = paste("grep -v '^#'", ref_gtf, "| grep -E '\t(gene|transcript|exon|CDS|cds|UTR|utr|gap)\t'"),
    sep = "\t", header = FALSE, quote = "", 
    col.names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  )
  
  fast_extract <- function(text, key) {
    pattern <- paste0(key, ' "([^"]*)"')
    res <- stri_match_first_regex(text, pattern)
    return(res[, 2])
  }
  
  gtf_df[, `:=`(
    gene_id = fast_extract(attributes, "gene_id"),
    gene_name = fast_extract(attributes, "gene_name"),
    transcript_id = fast_extract(attributes, "transcript_id"),
    gene_biotype = fast_extract(attributes, "gene_biotype")
  )]
  
  if (all(is.na(gtf_df$gene_biotype))) {
    gtf_df[, gene_biotype := fast_extract(attributes, "gene_type")]
  }
  data.table::setDF(gtf_df)
}

# 逻辑补全
gtf_df$gene_biotype[is.na(gtf_df$gene_biotype)] <- "protein_coding"
gtf_df$gene_name[is.na(gtf_df$gene_name)] <- gtf_df$gene_id[is.na(gtf_df$gene_name)]
gtf_df$tx_id <- gtf_df$transcript_id
gtf_df$tx_id[is.na(gtf_df$tx_id)] <- gtf_df$gene_id[is.na(gtf_df$tx_id)]

# 转换为 GRanges
gtf_filter <- makeGRangesFromDataFrame(gtf_df, keep.extra.columns = TRUE)
gtf_filter$type <- gtf_df$feature
# 清洗下划线
if (any(grepl("_", gtf_filter$gene_name, fixed = TRUE))) {
  gtf_filter$gene_name <- gsub("_", "-", gtf_filter$gene_name, fixed = TRUE)
}


obj[["ATAC"]] <- CreateChromatinAssay(counts = atac_data, sep = c(":", "-"), fragments = fragpath, annotation = gtf_filter)
obj


obj@meta.data$orig.ident <- rawname


# RNA cluster
DefaultAssay(obj) <- "RNA"
obj[['percent.mito']] <- PercentageFeatureSet(object = obj, pattern = '^(MT|mt|Mt)-')
obj@meta.data <- obj@meta.data %>% dplyr::rename("mito" = "percent.mito")
# RNA normalized
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# 健壮性改进：检查 mito 方差，若无方差则跳过回归
mito_values <- obj@meta.data$mito
mito_variance <- var(mito_values, na.rm = TRUE)
if (is.na(mito_variance) || mito_variance == 0) {
  cat("Warning: mito variance is zero or NA, skipping regression in ScaleData\n")
  obj <- ScaleData(obj)
} else {
  obj <- ScaleData(obj, vars.to.regress = "mito")
}

# 健壮性改进：动态调整 PCA 维数和邻居搜索维数
num_cells <- ncol(obj)
npcs_to_use <- min(50, num_cells - 1)
obj <- RunPCA(obj, npcs = npcs_to_use)

# RNA dimensionality reduction & clustering
dims_to_use <- 1:min(30, npcs_to_use)
# 健壮性改进：动态调整邻居数和聚类分辨率，防止低细胞量时产生过多碎片化 cluster
k_param <- min(20, num_cells - 1)
res_to_use <- ifelse(num_cells < 100, 0.2, 0.8)
obj <- FindNeighbors(obj, dims = dims_to_use, k.param = k_param)
obj <- FindClusters(obj, resolution = res_to_use)

# 动态计算 perplexity 和 neighbors 避免细胞数过少时报错
num_cells <- ncol(obj)
tsne_perplexity <- min(30, floor((num_cells - 1) / 3))
umap_neighbors <- min(30, num_cells - 1)
cat("RNA t-SNE perplexity set to:", tsne_perplexity, "\n")
cat("RNA UMAP n.neighbors set to:", umap_neighbors, "\n")

obj <- RunTSNE(obj, dims = dims_to_use, check_duplicates = FALSE, perplexity = tsne_perplexity)
obj <- RunUMAP(obj, dims = dims_to_use, n.neighbors = umap_neighbors)
# RNA tsne
tsne_loci <- as.data.frame(Embeddings(obj, reduction='tsne'))
tsne_loci <- cbind(tsne_loci, obj[[]])
write.table(tsne_loci, file='gex_tsne_umi.xls', 
            row.names=TRUE, 
            col.names=TRUE, 
            sep="\t", 
            quote=FALSE)


# diff table
features_df <- read.table(file.path(gex_matrix, 'features.tsv.gz'), sep="\t")
names(features_df)[1:2] <- c('Ensembl', 'gene')

# 健壮性改进：检查分群数量并安全处理 FindAllMarkers
if (length(unique(Idents(obj))) > 1) {
    obj.markers <- tryCatch({
        FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)
    }, error = function(e) {
        warning("FindAllMarkers failed: ", e$message)
        return(NULL)
    })
    
    if (!is.null(obj.markers) && nrow(obj.markers) > 0) {
        obj.markers <- obj.markers %>%
            left_join(features_df, by='gene') %>% 
            relocate(Ensembl, gene)
        obj.markers$Ensembl[is.na(obj.markers$Ensembl)] <- "na"
    } else {
        cat("No marker genes found.\n")
        obj.markers <- data.frame(Ensembl=character(), gene=character(), p_val=numeric(), avg_log2FC=numeric(), pct.1=numeric(), pct.2=numeric(), p_val_adj=numeric(), cluster=character())
    }
} else {
    cat("Only one cluster found, skipping FindAllMarkers.\n")
    obj.markers <- data.frame(Ensembl=character(), gene=character(), p_val=numeric(), avg_log2FC=numeric(), pct.1=numeric(), pct.2=numeric(), p_val_adj=numeric(), cluster=character())
}

write.table(obj.markers, file='gex_FindAllMarkers.xls', row.names=FALSE, sep="\t", quote=FALSE)


# ATAC cluster
# Assay ATAC count NS, TSS
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)
# ATAC normalized, dimensionality reduction & clustering
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# 动态调整 ATAC 维数
num_cells <- ncol(obj)
# LSI 第一维通常是测序深度，一般从第 2 维开始
atac_dims <- 2:min(30, num_cells - 1)
if (length(atac_dims) < 2) atac_dims <- 1:min(30, num_cells - 1)

umap_neighbors <- min(30, num_cells - 1)
obj <- RunUMAP(object = obj, reduction = 'lsi', dims = atac_dims, reduction.name = "atacumap", reduction.key = "atacumap_", n.neighbors = umap_neighbors)

# 动态计算 ATAC t-SNE perplexity
tsne_perplexity <- min(30, floor((num_cells - 1) / 3))
cat("ATAC t-SNE perplexity set to:", tsne_perplexity, "\n")
cat("ATAC UMAP n.neighbors set to:", umap_neighbors, "\n")

obj <- RunTSNE(obj, reduction = 'lsi', dims = atac_dims, reduction.name = "atactsne", reduction.key = "atactsne_", check_duplicates = FALSE, perplexity = tsne_perplexity)

# ATAC 也需要动态邻居数和分辨率
k_param_atac <- min(20, num_cells - 1)
res_to_use_atac <- ifelse(num_cells < 100, 0.2, 0.8)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = atac_dims, k.param = k_param_atac)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3, resolution = res_to_use_atac)
# ATAC tsne
tsne_loci <- as.data.frame(Embeddings(obj, reduction='atactsne'))
tsne_loci <- cbind(tsne_loci, obj[[]])
write.table(tsne_loci, file='atac_tsne_umi.xls', 
            row.names=TRUE, 
            col.names=TRUE, 
            sep="\t", 
            quote=FALSE)

saveRDS(obj,file=paste0(samplename,'.rds'))

# Linking peaks to genes
fa_data <- open(FaFile(ref_fa))
DefaultAssay(obj) <- "ATAC"

tryCatch({
  obj <- RegionStats(obj, genome = fa_data)
  obj <- LinkPeaks(
    object = obj, 
    peak.assay = "ATAC", 
    expression.assay = "RNA", 
    score_cutoff = 0.01, 
    distance = 1e6
  )
  saveRDS(obj,file='joint_peak_link_gene.rds')
  # count link
  linked_peaks <- Links(obj)
  total_links <- length(linked_peaks)
  # count link to gene
  linked_genes <- unique(linked_peaks$gene)
  total_linked_genes <- length(linked_genes)
  # count link to peaks
  linked_peaks_names <- unique(linked_peaks$peak)
  total_linked_peaks <- length(linked_peaks_names)
  # output result
  write.table(linked_peaks, file="linked_feature.xls", sep = "\t", quote = FALSE, row.names = FALSE)
}, error = function(e) {
  cat("The reference genome may be in NCBI format and can not calculate links. Please check!\n")
})
