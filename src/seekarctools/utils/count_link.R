suppressPackageStartupMessages({
  library(future)
  library(Signac)
  library(Seurat)
  library(dplyr)
  library(argparse)
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
parser$add_argument("--organism", help="human or mouse")
parser$add_argument("--anno_rds", help="Anno_EnsDb.rds")
args <- parser$parse_args()

organism=args$organism
if (organism == "human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
} else if (organism == "mouse") {
    library(BSgenome.Mmusculus.UCSC.mm10)
} else {
    stop("Warning: Not human or mouse, pipeline does not run astep4")
}

outdir=args$outdir
gex_matrix=args$gex_matrix
atac_matrix=args$atac_matrix
fragpath=args$fragpath
rawname=args$rawname
samplename=args$samplename
core=args$core
memory=args$memory
anno_rds=args$anno_rds

cat("available memory:", memory, "\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# multicore
options(future.globals.maxSize = as.integer(memory) * 1024 ^ 3)
plan("multicore", workers = as.integer(core))
plan()

gex_data <- Read10X(data.dir = gex_matrix)
gexobj <- CreateSeuratObject(counts = gex_data, assay = "RNA")
cat("------------gex------------------------\n")
gexobj

atac_data <- Read10X(data.dir = atac_matrix)
atacobj <- CreateSeuratObject(counts = atac_data, assay = "ATAC")
cat("------------atac------------------------\n")
atacobj

cat("------------joint------------------------\n")
# jointcb <- colnames(atacobj)
# obj <- subset(gexobj, cells = jointcb)
obj <- gexobj

annotation <- readRDS(anno_rds)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

obj[["ATAC"]] <- CreateChromatinAssay(counts = atac_data, sep = c(":", "-"), fragments = fragpath, annotation = annotation)
obj
cat("------------joint end------------------------\n")


cat("------------filter start------------------------\n")
# filter
DefaultAssay(obj) <- "ATAC"
features.keep <- as.character(seqnames(granges(obj))) %in% standardChromosomes(granges(obj))
obj.filter <- obj[features.keep, ]
obj[["ATAC"]] <- obj.filter[["ATAC"]]
# obj[["ATAC"]] <- subset(obj[["ATAC"]], features = rownames(obj[["ATAC"]])[features.keep]) # seurat 5 
obj@meta.data$orig.ident <- rawname
cat("------------filter end------------------------\n")


cat("------------RNA cluster start------------------------\n")
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
cat("------------RNA cluster end------------------------\n")


# diff table
features_df <- read.table(file.path(gex_matrix, 'features.tsv.gz'), sep="\t")
names(features_df)[1:2] <- c('Ensembl', 'gene')
obj.markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE) %>%
    left_join(features_df, by='gene') %>% relocate(Ensembl, gene)
obj.markers$Ensembl[is.na(obj.markers$Ensembl)] <- "na"
write.table(obj.markers, file='gex_FindAllMarkers.xls', row.names=FALSE, sep="\t", quote=FALSE)


cat("------------ATAC cluster start------------------------\n")
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
cat("------------ATAC cluster end------------------------\n")

cat("------------WNN umap start------------------------\n")
## WNN
DefaultAssay(obj) <- "RNA"
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnnumap",
  reduction.key = "wnnumap_",
  assay = "RNA",
  verbose = TRUE
)
cat("------------WNN tsne start------------------------\n")
obj <- RunTSNE(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnntsne",
  reduction.key = "wnntsne_",
  assay = "RNA",
  verbose = TRUE,
  check_duplicates = FALSE
)
cat("------------WNN cluster end------------------------\n")
saveRDS(obj,file=paste0(samplename,'.rds'))

cat("------------Linking peaks to genes start------------------------\n")
# Linking peaks to genes
DefaultAssay(obj) <- "ATAC"
if (organism == "human") {
    obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
} else {
    obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10)
}
obj <- LinkPeaks(object = obj,peak.assay = "ATAC",expression.assay = "RNA")

saveRDS(obj,file='joint_peak_link_gene.rds')
cat("------------Linking peaks to genes end------------------------\n")

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
cat("total links:", total_links, "\n")
cat("total linked genes:", total_linked_genes, "\n")
cat("total linked peaks:", total_linked_peaks, "\n")

write.table(linked_peaks, file="linked_feature.xls", sep = "\t", quote = FALSE, row.names = FALSE)


