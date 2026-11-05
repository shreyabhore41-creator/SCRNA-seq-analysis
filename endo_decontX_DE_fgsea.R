## Re-run library-level DE after removing ambient RNA, then pathway analysis
## Run from your project folder (same one as your main script).
## Needs obj_annotated.rds (your annotated Seurat object) and the clean metadata.

# BiocManager::install(c("celda", "fgsea", "msigdbr"))   # first time only
library(Seurat); library(SingleCellExperiment); library(celda)
library(DESeq2); library(dplyr); library(Matrix)
library(fgsea); library(msigdbr); library(ggplot2)
set.seed(42)

obj  <- readRDS("obj_annotated.rds")
meta <- read.csv("GSE203191_metadata_clean.csv", stringsAsFactors = FALSE)

## ---------------------------------------------------------------
## 1. Ambient RNA removal (DecontX), per library, using cell types
## ---------------------------------------------------------------
obj <- JoinLayers(obj)
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce <- decontX(sce, z = sce$celltype, batch = sce$library)

# How much contamination per library? (useful number for your report)
print(aggregate(decontX_contamination ~ library, data = as.data.frame(colData(sce)), mean))

clean_counts <- round(decontXcounts(sce))
colnames(clean_counts) <- colnames(obj)
obj[["decont"]] <- CreateAssayObject(counts = clean_counts)

## Genes to exclude from DE: haemoglobin + mitochondrial
bad_genes <- grep("^HB[ABDEGMQZ]|^HB[AB][0-9]?$|^MT-", rownames(obj), value = TRUE)
bad_genes <- union(bad_genes, c("HBB", "HBA1", "HBA2", "HBD"))

## ---------------------------------------------------------------
## 2. Library-level design table (same logic as your original script)
## ---------------------------------------------------------------
meta$library <- sub("_filtered.*", "", meta$processed_data_file)
frac <- as.data.frame.matrix(prop.table(table(meta$library, meta$phenotype), 1))
colnames(frac) <- make.names(colnames(frac))
ctrl <- grep("control", colnames(frac), ignore.case = TRUE, value = TRUE)
sym  <- grep("sympt",  colnames(frac), ignore.case = TRUE, value = TRUE)
case <- setdiff(colnames(frac), c(ctrl, sym))      # Diagnosed
print(list(ctrl = ctrl, sym = sym, case = case))   # each should be ONE name
frac$library <- rownames(frac)
lib_info <- frac[frac$library != "wholeME_MEtissue", ]   # unclear library left out
lib_info$prep <- ifelse(grepl("^wholeME[0-9]", lib_info$library), "whole_ME", "ME_tissue")

## ---------------------------------------------------------------
## 3. Pseudobulk DE on decontaminated counts
## ---------------------------------------------------------------
counts <- GetAssayData(obj, assay = "decont", layer = "counts")
counts <- counts[!rownames(counts) %in% bad_genes, ]

run_lib_de <- function(ct, effect, min_cells = 30) {
  md <- obj@meta.data
  n_per <- table(md$library[md$celltype == ct])
  ok <- intersect(names(n_per)[n_per >= min_cells], lib_info$library)
  if (length(ok) < 7) return(NULL)
  keep  <- md$celltype == ct & md$library %in% ok
  cells <- rownames(md)[keep]
  lib   <- factor(md$library[keep])
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
run_safe <- function(ct, effect) tryCatch(run_lib_de(ct, effect), error = function(e) NULL)

cts <- names(which(table(obj$celltype) >= 500))
cts <- setdiff(cts, "Mixed/unclear")
res_case <- bind_rows(lapply(cts, run_safe, effect = case))
res_sym  <- bind_rows(lapply(cts, run_safe, effect = sym))
write.csv(res_case, "results/libDE_decont_diagnosed_vs_control.csv", row.names = FALSE)
write.csv(res_sym,  "results/libDE_decont_symptomatic_vs_control.csv", row.names = FALSE)

# Did the contamination genes go away?
res_case %>% filter(padj < 0.1) %>% count(celltype)
res_case %>% filter(gene %in% c("LUM", "DCN", "RBP1")) %>%
  select(celltype, gene, log2FoldChange, pvalue, padj)

## ---------------------------------------------------------------
## 4. fgsea on the Wald statistic (per cell type, Hallmark sets)
## ---------------------------------------------------------------
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

run_gsea <- function(df) {
  df <- df %>% filter(!is.na(stat)) %>% group_by(gene) %>% slice_max(abs(stat), n = 1) %>% ungroup()
  ranks <- setNames(df$stat, df$gene)
  fgsea(pathways, ranks, minSize = 15, maxSize = 500) %>% arrange(padj)
}

gsea_all <- bind_rows(lapply(split(res_case, res_case$celltype), function(d) {
  r <- run_gsea(d); r$celltype <- d$celltype[1]; r
}))
gsea_all$leadingEdge <- sapply(gsea_all$leadingEdge, paste, collapse = ";")
write.csv(gsea_all, "results/fgsea_hallmark_diagnosed_vs_control.csv", row.names = FALSE)
print(gsea_all %>% filter(padj < 0.1) %>% select(celltype, pathway, NES, padj), n = 50)

## ---------------------------------------------------------------
## 5. Volcano plots (labels use raw p because power is low)
## ---------------------------------------------------------------
v <- res_case %>% filter(!is.na(pvalue)) %>%
  mutate(sig = ifelse(!is.na(padj) & padj < 0.1, "padj<0.1", "ns"))
pv <- ggplot(v, aes(log2FoldChange, -log10(pvalue), colour = sig)) +
  geom_point(size = 0.6, alpha = 0.6) + facet_wrap(~celltype) + theme_bw() +
  scale_colour_manual(values = c(ns = "grey70", "padj<0.1" = "firebrick")) +
  labs(x = "log2FC per unit fraction of Diagnosed in library (NOT a standard log2FC)")
ggsave("figures/volcano_decont_diagnosed.png", pv, width = 12, height = 9, dpi = 300)
