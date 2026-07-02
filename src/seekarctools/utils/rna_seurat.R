suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(argparse)
  library(future)
})

parser <- ArgumentParser()
parser$add_argument("--gex_matrix", required = TRUE, help = "RNA filtered_feature_bc_matrix")
parser$add_argument("--samplename", required = TRUE, help = "sample name")
parser$add_argument("--outdir", required = TRUE, help = "output directory")
parser$add_argument("--core", default = 4, help = "parallel cores")
parser$add_argument("--memory", default = 32, help = "memory GB")
args <- parser$parse_args()

gex_matrix <- args$gex_matrix
samplename <- args$samplename
outdir <- args$outdir
core <- as.integer(args$core)
memory <- as.integer(args$memory)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

options(future.globals.maxSize = memory * 1024 ^ 3)
options(future.rng.onMisuse = "ignore")
plan("multicore", workers = core)

cat("[RNA RDS] read matrix:", gex_matrix, "\n")
gex_data <- Read10X(data.dir = gex_matrix)

if (is.list(gex_data)) {
  if ("Gene Expression" %in% names(gex_data)) {
    gex_data <- gex_data[["Gene Expression"]]
  } else {
    gex_data <- gex_data[[1]]
  }
}

obj <- CreateSeuratObject(counts = gex_data, assay = "RNA", project = samplename)
obj@meta.data$orig.ident <- samplename

num_cells <- ncol(obj)
num_features <- nrow(obj)
cat("[RNA RDS] cells:", num_cells, "\n")
cat("[RNA RDS] features:", num_features, "\n")

obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^(MT|mt|Mt)-")
obj@meta.data <- obj@meta.data %>% dplyr::rename("mito" = "percent.mito")

if (num_cells >= 3 && num_features >= 50) {
  cat("[RNA RDS] normalize and cluster\n")

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000, num_features))

  mito_values <- obj@meta.data$mito
  mito_variance <- var(mito_values, na.rm = TRUE)

  if (is.na(mito_variance) || mito_variance == 0) {
    cat("[RNA RDS] mito variance is zero or NA, skip regression\n")
    obj <- ScaleData(obj)
  } else {
    obj <- ScaleData(obj, vars.to.regress = "mito")
  }

  npcs_to_use <- min(50, num_cells - 1, num_features - 1)
  obj <- RunPCA(obj, npcs = npcs_to_use)

  dims_to_use <- 1:min(30, npcs_to_use)
  k_param <- min(20, num_cells - 1)
  res_to_use <- ifelse(num_cells < 100, 0.2, 0.8)

  obj <- FindNeighbors(obj, dims = dims_to_use, k.param = k_param)
  obj <- FindClusters(obj, resolution = res_to_use)

  if (num_cells >= 4) {
    tsne_perplexity <- min(30, floor((num_cells - 1) / 3))
    if (tsne_perplexity >= 1) {
      cat("[RNA RDS] tSNE perplexity:", tsne_perplexity, "\n")
      obj <- RunTSNE(obj, dims = dims_to_use, check_duplicates = FALSE, perplexity = tsne_perplexity)

      tsne_loci <- as.data.frame(Embeddings(obj, reduction = "tsne"))
      tsne_loci <- cbind(tsne_loci, obj[[]])
      write.table(
        tsne_loci,
        file = "gex_tsne_umi.xls",
        row.names = TRUE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
      )
    }
  }

  if (num_cells >= 3) {
    umap_neighbors <- min(30, num_cells - 1)
    cat("[RNA RDS] UMAP n.neighbors:", umap_neighbors, "\n")
    obj <- RunUMAP(obj, dims = dims_to_use, n.neighbors = umap_neighbors)
  }

  if (length(unique(Idents(obj))) > 1) {
    obj.markers <- tryCatch({
      FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)
    }, error = function(e) {
      warning("FindAllMarkers failed: ", e$message)
      return(NULL)
    })

    if (is.null(obj.markers) || nrow(obj.markers) == 0) {
      obj.markers <- data.frame(
        gene = character(),
        p_val = numeric(),
        avg_log2FC = numeric(),
        pct.1 = numeric(),
        pct.2 = numeric(),
        p_val_adj = numeric(),
        cluster = character()
      )
    }
  } else {
    cat("[RNA RDS] only one cluster, skip FindAllMarkers\n")
    obj.markers <- data.frame(
      gene = character(),
      p_val = numeric(),
      avg_log2FC = numeric(),
      pct.1 = numeric(),
      pct.2 = numeric(),
      p_val_adj = numeric(),
      cluster = character()
    )
  }

  write.table(
    obj.markers,
    file = "gex_FindAllMarkers.xls",
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
} else {
  cat("[RNA RDS] too few cells or features, only save raw Seurat object\n")
}

rds_file <- paste0(samplename, ".rds")
saveRDS(obj, file = rds_file)
cat("[RNA RDS] saved:", file.path(outdir, rds_file), "\n")