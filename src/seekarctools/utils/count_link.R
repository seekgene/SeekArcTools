library(future)
library(Signac)
library(Seurat)
library(dplyr)
library(BiocParallel)
library(argparse)
parser = ArgumentParser()
parser$add_argument("--gex_matrix", help="gex. step3/filtered_feature_bc_matrix")
parser$add_argument("--atac_matrix", help="atac. step3/filter_peaks_bc_matrix")
parser$add_argument("--fragpath", help="atac. step3/asample_fragments.tsv.gz")
parser$add_argument("--outdir", help="outdir")
parser$add_argument("--core", help="Parallel running cores")
parser$add_argument("--species", help="human or mouse")
parser$add_argument("--anno_rds", help="Anno_EnsDb_Hsapiens_v86.rds or Anno_EnsDb_Mmusculus_v79.rds")
parser$add_argument("--memory", help="Memory usage")
args <- parser$parse_args()

species=args$species
if (species == "human") {
    library(BSgenome.Hsapiens.UCSC.hg38)
} else if (species == "mouse") {
    library(BSgenome.Mmusculus.UCSC.mm10)
} else {
    stop("Not human or mouse, please enter: 'human' or 'mouse'")
}

outdir=args$outdir
gex_matrix=args$gex_matrix
atac_matrix=args$atac_matrix
fragpath=args$fragpath
core=args$core
anno_rds=args$anno_rds
memory=args$memory

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

# jointcb <- colnames(atacobj)
# obj <- subset(gexobj, cells = jointcb)
cat("------------joint------------------------\n")
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
saveRDS(obj,file='filter_peaks_bc_matrix.rds')
cat("------------filter end------------------------\n")


cat("------------RNA cluster start------------------------\n")
# RNA normalized
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj)
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
obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30, reduction.name = "umapATAC", reduction.key = "umapATAC_")
obj <- RunTSNE(obj, reduction = 'lsi', dims = 2:30, reduction.name = "tsneATAC", reduction.key = "tsneATAC_")
obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
# ATAC tsne
tsne_loci <- as.data.frame(Embeddings(obj, reduction='tsneATAC'))
tsne_loci <- cbind(tsne_loci, obj[[]])
write.table(tsne_loci, file='atac_tsne_umi.xls', 
            row.names=TRUE, 
            col.names=TRUE, 
            sep="\t", 
            quote=FALSE)
cat("------------ATAC cluster end------------------------\n")

cat("------------WNN cluster start------------------------\n")
## WNN
DefaultAssay(obj) <- "SCT"
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
  reduction.name = "umapWNN",
  reduction.key = "umapWNN_",
  assay = "RNA",
  verbose = TRUE
)
cat("------------WNN cluster end------------------------\n")


cat("------------Linking peaks to genes start------------------------\n")
# Linking peaks to genes
DefaultAssay(obj) <- "ATAC"
if (species == "human") {
    obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
} else {
    obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10)
}
obj <- LinkPeaks(object = obj,peak.assay = "ATAC",expression.assay = "SCT")

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


