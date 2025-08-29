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
plan("multicore", workers = as.integer(core))
plan()

gex_data <- Read10X(data.dir = gex_matrix)
obj <- CreateSeuratObject(counts = gex_data, assay = "RNA")
atac_data <- Read10X(data.dir = atac_matrix)

# read gtf
gtf_data <- rtracklayer::import(ref_gtf, format = "gtf")
# make GRanges obj
mc <- mcols(gtf_data)

if ("gene_type" %in% colnames(mc)) {
      mcols(gtf_data)$gene_biotype <- mcols(gtf_data)$gene_type
    } else if ("biotype" %in% colnames(mc)) {
      mcols(gtf_data)$gene_biotype <- mcols(gtf_data)$biotype
    } else if ("gene_biotype" %in% colnames(mc)) {
      cat("There is a gene_biotype column in GTF.\n")
    } else {
      mcols(gtf_data)$gene_biotype <- rep("protein_coding", length(mc$gene_id))
    }

# Check and replace the NA value in gene_fiotype with 'proteino_comding'
na_indices <- is.na(mcols(gtf_data)$gene_biotype)
if (any(na_indices)) {
  mcols(gtf_data)$gene_biotype[na_indices] <- "protein_coding"
  cat("Replaced NA values in gene_biotype with 'protein_coding'.\n")
}

if (!("gene_name" %in% colnames(mc))) {
  mcols(gtf_data)$gene_name <- mcols(gtf_data)$gene_id
  cat("gene_name do not exist, use gene_id instead\n")
}

# Check and replace the NA value in gene_name with gene_id
missing_names <- is.na(mcols(gtf_data)$gene_name)
if (any(missing_names)) {
  mcols(gtf_data)$gene_name[missing_names] <- mcols(gtf_data)$gene_id[missing_names]
  cat("gene_name exist NA value, use gene_id instead\n")
}

if ("tx_id" %in% colnames(mc)) {
      cat("There is a tx_id column in GTF.\n")
    } else if ("transcript_id" %in% colnames(mc)) {
      mcols(gtf_data)$tx_id <- mcols(gtf_data)$transcript_id
    } else {
      mcols(gtf_data)$tx_id <- mcols(gtf_data)$gene_id
      cat("transcript_id do not exist, use gene_id instead\n")
    }

keep_cols <- c("tx_id", "gene_name", "gene_id", "gene_biotype", "type")
mcols(gtf_data) <- mcols(gtf_data)[, keep_cols]
gtf_filter <- gtf_data[gtf_data$type %in% c("CDS", "cds", "UTR", "utr", "exon", "gap")]


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
obj <- ScaleData(obj, vars.to.regress = "mito")
obj <- RunPCA(obj)
# RNA dimensionality reduction & clustering
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8)
obj <- RunTSNE(obj, dims = 1:30,check_duplicates = FALSE)
obj <- RunUMAP(obj, dims = 1:30)
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
obj.markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE) %>%
    left_join(features_df, by='gene') %>% relocate(Ensembl, gene)
obj.markers$Ensembl[is.na(obj.markers$Ensembl)] <- "na"
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
obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30, reduction.name = "atacumap", reduction.key = "atacumap_")
obj <- RunTSNE(obj, reduction = 'lsi', dims = 2:30, reduction.name = "atactsne", reduction.key = "atactsne_", check_duplicates = FALSE)
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
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
  obj <- LinkPeaks(object = obj, peak.assay = "ATAC", expression.assay = "RNA")
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
