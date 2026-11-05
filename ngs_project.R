###############################################################################
# GSE203191 ENDOMETRIOSIS scRNA-seq ANALYSIS
# Cleaned copy of the user's original script.
#
# Only exploratory/debugging commands, repeated installation/loading steps,
# duplicate analysis blocks, and patient-level analyses not supported by the
# pooled-library design have been removed. Retained analysis code and
# thresholds are otherwise kept from the original script.
#
# DecontX + corrected library-level DE + fgsea are maintained separately in
# endo_decontX_DE_fgsea.R (File 1), so that block is not duplicated here.
###############################################################################

library(Seurat); library(harmony); library(dplyr); library(ggplot2)
setwd("/home/student/Downloads/NGS_mini project")

library(harmony)
library(dplyr)
library(ggplot2)
library(Matrix); library(patchwork); library(SingleCellExperiment)
library(scDblFinder); library(DESeq2); library(speckle)
set.seed(42)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

meta <- utils::read.csv(
  "GSE203191_metadata_clean.csv",
  stringsAsFactors = FALSE
)

data_dir <- "/home/student/Downloads/NGS_mini project/GSE203191_RAW"
files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)
length(files)
basename(files)

# Correct data directory
data_dir <- "/home/student/Downloads/NGS_mini project/GSE203191_RAW"

# Get full paths of all H5 files
files <- list.files(
  data_dir,
  pattern = "\\.h5$",
  full.names = TRUE
)

# Check files exist
stopifnot(length(files) == 11)
stopifnot(all(file.exists(files)))

# Load expression matrices into Seurat
objs <- lapply(files, function(f) {
  
  m <- Seurat::Read10X_h5(f)
  
  obj <- Seurat::CreateSeuratObject(
    counts = m,
    project = basename(f),
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample_file <- basename(f)
  
  return(obj)
})


for (i in seq_along(objs)) {
  file_name <- basename(files[i])
  file_key <- sub("^GSM[0-9]+_", "", file_name)
  
  rows <- which(meta$processed_data_file == file_key)
  
  print(file_name)
  print(meta[rows, c("subject_id", "phenotype", "tissue")])
}


head(meta[, c("title", "subject_id", "processed_data_file", "phenotype", "tissue")], 15)

table(meta$tissue, meta$phenotype)
sapply(objs, dim)

sapply(objs, function(x) inherits(x, "Seurat"))

obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))

dim(obj)

obj$library <- sub("_filtered.*", "", sub("^GSM[0-9]+_", "", obj$sample_file))
table(obj$library)

# 1. Preparation type (the combined-name file is unclear, so we flag it)
obj$prep <- ifelse(grepl("^wholeME[0-9]", obj$library), "whole_ME",
                   ifelse(grepl("^MEtissue", obj$library), "ME_tissue", "unclear"))
table(obj$library, obj$prep)

# 2. QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]")

# 3. Readable QC plot (pt.size = 0 removes the black dots)
p <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "library", pt.size = 0, ncol = 3)
p
# 4. Filter (this is my reading of the paper's thresholds: nUMI 500-50,000,
#    under 6,000 genes, under 10% mitochondrial reads)
obj <- subset(obj, subset = nCount_RNA > 500 & nCount_RNA < 50000 &
                nFeature_RNA < 6000 & percent.mt < 10)
table(obj$library)   # cells left per library


# 1. Doublets, run separately per library
library(scDblFinder); library(SingleCellExperiment)
obj <- JoinLayers(obj)          # if this errors, tell me and skip it
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce, samples = "library")
obj$doublet_class <- sce$scDblFinder.class
table(obj$library, obj$doublet_class)
obj <- subset(obj, subset = doublet_class == "singlet")


# 2. Normalize and reduce
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)

dim(obj)                    # cells after removing doublets
table(obj$doublet_class)    # should show only "singlet"

# 1. Remove library-to-library differences
library(harmony)
obj <- RunHarmony(obj, group.by.vars = "library")

# 2. UMAP and clustering using 20 dimensions
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)

# 3. Look at it and save
u <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE) |
  DimPlot(obj, group.by = "library")
u
ggsave("figures/umap.png", u, width = 14, height = 6, dpi = 300)
saveRDS(obj, "obj_after_umap.rds")

DimPlot(obj, group.by = "library", shuffle = TRUE)
table(obj$seurat_clusters, obj$library)

markers <- c("COL1A1","DCN","PDGFRB",       # stromal
             "EPCAM","KRT8","PAEP",          # epithelial
             "PECAM1","VWF",                 # endothelial
             "RGS5","ACTA2",                 # perivascular/mural
             "NKG7","GNLY","KLRC1",          # NK
             "CD3E","IL7R",                  # T cells
             "MS4A1","CD79A",                # B cells
             "CD14","LYZ","CD68")            # myeloid
d <- DotPlot(obj, features = markers, group.by = "seurat_clusters") + RotatedAxis()
d
ggsave("figures/marker_dotplot.png", d, width = 14, height = 6, dpi = 300)

all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all_markers, "results/cluster_markers.csv", row.names = FALSE)


# Fill this in by hand after reading the dotplot and marker table
cluster_labels <- c(
  "0" = "Stromal"
  # "1" = "...", "2" = "...", ...
)
obj$celltype <- unname(cluster_labels[as.character(obj$seurat_clusters)])

## 6. Map onto a tissue reference (optional but central to your question)
# Needs a prepared reference Seurat object (e.g. from Tan 2022 or HECA) with
# normalised data, PCA, and a cell-state label column.
# ref <- readRDS("reference/endometrium_reference.rds")
# anchors <- FindTransferAnchors(reference = ref, query = obj,
#                                dims = 1:30, reference.reduction = "pca")
# pred <- TransferData(anchorset = anchors, refdata = ref$celltype, dims = 1:30)
# obj <- AddMetaData(obj, pred)   # adds predicted.id and prediction.score.max
# table(obj$predicted.id, obj$group)

rm(cluster_labels)
top5 <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
print(top5, n = 90)
saveRDS(obj, "obj_clustered.rds")

# 1. Check the marker results exist
dim(all_markers)        # should have thousands of rows
head(all_markers)

# 2. One compact line of top genes per cluster (easy to paste or screenshot)
library(dplyr)
top <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 8) %>%
  summarise(genes = paste(gene, collapse = ", "))
print(as.data.frame(top), right = FALSE)
# 3. Save it
write.csv(top, "results/top_markers_per_cluster.csv", row.names = FALSE)

# 1. Proportion of each cell type in each library
comp <- obj@meta.data %>%
  count(library, prep, celltype) %>%
  group_by(library) %>%
  mutate(prop = n / sum(n)) %>% ungroup()
write.csv(comp, "results/celltype_by_library.csv", row.names = FALSE)

comp <- obj@meta.data %>%
  dplyr::count(library, prep, celltype) %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()
write.csv(comp, "results/celltype_by_library.csv", row.names = FALSE)

FeaturePlot(obj, c("COL1A1", "DCN", "HBB", "PTPRC"), ncol = 2)

cp <- ggplot(comp, aes(prep, prop, colour = prep)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
  facet_wrap(~celltype, scales = "free_y") + theme_bw()
ggsave("figures/composition_by_prep.png", cp, width = 12, height = 8, dpi = 300)

fp <- FeaturePlot(obj, c("COL1A1", "DCN", "HBB", "PTPRC"), ncol = 2)
ggsave("figures/featureplot_stromal_check.png", fp, width = 9, height = 8, dpi = 300)

obj$celltype[obj$seurat_clusters == 8] <- "Mixed/unclear"
table(obj$celltype)
saveRDS(obj, "obj_annotated.rds")
l <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE)
ggsave("figures/umap_celltypes.png", l, width = 9, height = 7, dpi = 300)

FeaturePlot(obj, c("CD3E", "NKG7", "CD14", "MS4A1"), ncol = 2)

VlnPlot(obj, c("nFeature_RNA", "percent.mt", "PTPRC", "COL1A1", "LUM", "HBB"),
        group.by = "seurat_clusters", pt.size = 0, ncol = 3)
obj$celltype[obj$seurat_clusters == 8] <- "Mixed/unclear"
saveRDS(obj, "obj_annotated.rds")

sub <- obj@meta.data %>% dplyr::filter(prep %in% c("whole_ME", "ME_tissue"))
prop_res <- propeller(clusters = sub$celltype, sample = sub$library, group = sub$prep)
prop_res
write.csv(prop_res, "results/propeller_wholeME_vs_tissue.csv")

meta$library <- sub("_filtered.*", "", meta$processed_data_file)
table(meta$library, meta$phenotype)
cp <- as.data.frame(prop.table(table(library = obj$library, celltype = obj$celltype), 1),
                    responseName = "prop")


frac <- as.data.frame.matrix(prop.table(table(meta$library, meta$phenotype), 1))
colnames(frac) <- make.names(colnames(frac))
frac$library <- rownames(frac)

ctrl <- grep("control", names(frac), ignore.case = TRUE, value = TRUE)
sym  <- grep("sympt",  names(frac), ignore.case = TRUE, value = TRUE)
case <- setdiff(names(frac), c(ctrl, sym, "library"))
print(c(ctrl = ctrl, sym = sym, case = case))   # check these are right

dat <- merge(cp, frac, by = "library")
dat <- dat[dat$library != "wholeME_MEtissue", ]       # unclear library
dat$prep <- ifelse(grepl("^wholeME", dat$library), "whole_ME", "ME_tissue")


res <- do.call(rbind, lapply(split(dat, dat$celltype), function(d) {
  m <- lm(reformulate(c(case, sym, "prep"), "prop"), data = d)
  s <- summary(m)$coefficients
  data.frame(celltype = d$celltype[1],
             case_vs_control = s[case, "Estimate"], p_case = s[case, "Pr(>|t|)"],
             sym_vs_control  = s[sym,  "Estimate"], p_sym  = s[sym,  "Pr(>|t|)"])
}))
res$FDR_case <- p.adjust(res$p_case, "BH")
res[order(res$p_case), ]
write.csv(res, "results/library_level_case_vs_control.csv", row.names = FALSE)

######################################################################################
meta$library <- sub("_filtered.*", "", meta$processed_data_file)
frac <- as.data.frame.matrix(prop.table(table(meta$library, meta$phenotype), 1))
colnames(frac) <- make.names(colnames(frac))
ctrl <- grep("control", colnames(frac), ignore.case = TRUE, value = TRUE)
sym  <- grep("sympt",  colnames(frac), ignore.case = TRUE, value = TRUE)
case <- setdiff(colnames(frac), c(ctrl, sym))
print(list(ctrl = ctrl, sym = sym, case = case))   # each should be ONE name
frac$library <- rownames(frac)
lib_info <- frac[frac$library != "wholeME_MEtissue", ]   # unclear library left out
lib_info$prep <- ifelse(grepl("^wholeME[0-9]", lib_info$library), "whole_ME", "ME_tissue")
print(lib_info)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

run_lib_de <- function(ct, effect, min_cells = 30) {
  md <- obj@meta.data
  n_per <- table(md$library[md$celltype == ct])
  ok <- intersect(names(n_per)[n_per >= min_cells], lib_info$library)
  if (length(ok) < 7) return(NULL)
  keep <- md$celltype == ct & md$library %in% ok
  cells <- rownames(md)[keep]
  lib <- factor(md$library[keep])
  pb <- counts[, cells] %*% Matrix::sparse.model.matrix(~ 0 + lib)
  colnames(pb) <- sub("^lib", "", colnames(pb))
  cd <- lib_info[match(colnames(pb), lib_info$library), ]
  rownames(cd) <- cd$library
  cd$prep <- factor(cd$prep)
  des <- reformulate(c(if (nlevels(cd$prep) > 1) "prep", case, sym))
  mm <- model.matrix(des, cd)
  if (qr(mm)$rank < ncol(mm)) return(NULL)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(pb)), cd, design = des)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds, quiet = TRUE)
  out <- as.data.frame(results(dds, name = effect))
  out$gene <- rownames(out); out$celltype <- ct
  out$effect <- effect; out$n_libs <- length(ok)
  out
}


library(DESeq2); library(dplyr)
cts <- names(which(table(obj$celltype) >= 500))
cts <- setdiff(cts, "Mixed/unclear")
run_safe <- function(ct, effect) tryCatch(run_lib_de(ct, effect), error = function(e) NULL)
res_case <- bind_rows(lapply(cts, run_safe, effect = case))
res_sym  <- bind_rows(lapply(cts, run_safe, effect = sym))
table(res_case$celltype)      # which cell types had enough libraries
write.csv(res_case, "results/libDE_case_vs_control.csv", row.names = FALSE)
write.csv(res_sym,  "results/libDE_symptomatic_vs_control.csv", row.names = FALSE)
top <- res_case %>% filter(!is.na(pvalue)) %>%
  group_by(celltype) %>% slice_min(pvalue, n = 10) %>% ungroup()
print(as.data.frame(top[, c("celltype","gene","log2FoldChange","pvalue","padj")]))

chk <- c("IGFBP1","LEFTY2","IL11","MGP","FN1","MMP1")
res_case %>% filter(celltype %in% c("Stromal"), gene %in% chk) %>%
  select(gene, log2FoldChange, pvalue)
