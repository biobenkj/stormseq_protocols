## principle curves with slingshot
suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  library(tximeta)
  library(rjson)
  library(SingleCellExperiment)
  library(scater)
  library(slingshot)
})

# load data
salmon_quants_hg38 <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")
salmon_quants_hg38.hiseq <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")

# collapse clusters
library(dplyr)
salmon_quants_hg38$cluster.collapse <- salmon_quants_hg38$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "Ciliated Intermediate",
    . == 2 ~ "Secretory Intermediate",
    . == 3 ~ "UCFP",
    . == 4 ~ "Ciliated",
    . == 5 ~ "Branch",
    . == 6 ~ "Immune",
    . == 7 ~ "Secretory",
    . == 8 ~ "Secretory"
  ))
salmon_quants_hg38$cluster.collapse <- salmon_quants_hg38$cluster.collapse$cluster_cat

salmon_quants_hg38.hiseq$cluster.collapse <- salmon_quants_hg38.hiseq$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "UCFP",
    . == 2 ~ "Secretory",
    . == 3 ~ "Ciliated",
    . == 4 ~ "Secretory",
    . == 5 ~ "Secretory",
    . == 6 ~ "Ciliated"
  ))
salmon_quants_hg38.hiseq$cluster.collapse <- salmon_quants_hg38.hiseq$cluster.collapse$cluster_cat


# remove immune for now to match densmap embeddings/trajectory inference
salmon_quants_hg38.noimmune <- salmon_quants_hg38[,!salmon_quants_hg38$cluster.collapse %in% "Immune"]

# filter out anything with CD45 expression :)
plotReducedDim(salmon_quants_hg38.hiseq, dimred = "PCAsub", colour_by = "cd45")
salmon_quants_hg38.hiseq.noimmune <- salmon_quants_hg38.hiseq[,!salmon_quants_hg38.hiseq$cd45 > 1.0]


# quick plot to make sure things are copacetic 
plotReducedDim(salmon_quants_hg38.noimmune, dimred = "PCAsub", colour_by = "cluster.collapse")
plotReducedDim(salmon_quants_hg38.hiseq.noimmune, dimred = "PCAsub", colour_by = "cluster.collapse")

# principle curve fitting
set.seed(42)
library(RColorBrewer)
salmon_quants_hg38.noimmune.lin <- getLineages(salmon_quants_hg38.noimmune,
                                               clusterLabels = "cluster.collapse",
                                               reducedDim = "PCAsub",
                                               start.clus = "UCFP")
salmon_quants_hg38.noimmune.lin.curves <- getCurves(salmon_quants_hg38.noimmune.lin)
plot(reducedDims(salmon_quants_hg38.noimmune)$PCAsub, col = brewer.pal(9,"Set1")[factor(salmon_quants_hg38.noimmune$cluster.collapse)], asp = 1, pch = 16)
lines(SlingshotDataSet(salmon_quants_hg38.noimmune.lin.curves), lwd = 3, col = 'black')


curves.nova <- slingCurves(SlingshotDataSet(salmon_quants_hg38.noimmune.lin.curves))

salmon_quants_hg38.hiseq.lin <- getLineages(salmon_quants_hg38.hiseq.noimmune,
                                               clusterLabels = "cluster.collapse",
                                               reducedDim = "PCAsub",
                                               start.clus = "UCFP")
salmon_quants_hg38.hiseq.lin.curves <- getCurves(salmon_quants_hg38.hiseq.lin)
plot(reducedDims(salmon_quants_hg38.hiseq.noimmune)$PCAsub, col = brewer.pal(9,"Set1")[factor(salmon_quants_hg38.hiseq.noimmune$cluster.collapse)], asp = 1, pch = 16)
lines(SlingshotDataSet(salmon_quants_hg38.hiseq.lin.curves), lwd = 3, col = 'black')
curves.hiseq <- slingCurves(SlingshotDataSet(salmon_quants_hg38.hiseq.lin.curves))


# extract curves for plotting
curve_df.nova <- lapply(seq_along(curves.nova), function(i) {
  curve_i <- curves.nova[[i]]
  #browser()
  data.frame(
    Dim1 = curve_i$s[,1],
    Dim2 = curve_i$s[,2],
    #pseudotime = curve_i$lambda,
    lineage = paste0("Lineage ", i)
  )
}) %>% bind_rows()

curve_df.hiseq <- lapply(seq_along(curves.hiseq), function(i) {
  curve_i <- curves.hiseq[[i]]
  #browser()
  data.frame(
    Dim1 = curve_i$s[,1],
    Dim2 = curve_i$s[,2],
    #pseudotime = curve_i$lambda,
    lineage = paste0("Lineage ", i)
  )
}) %>% bind_rows()

# get the pseudotime values
sling_time.nova <- slingAvgPseudotime(salmon_quants_hg38.noimmune.lin.curves)
sling_time.hiseq <- slingAvgPseudotime(salmon_quants_hg38.hiseq.lin.curves)

# rescale to 0-1
sling_time.nova <- scales::rescale(sling_time.nova, to = c(0,1))
sling_time.hiseq <- scales::rescale(sling_time.hiseq, to = c(0,1))

# updates to tack on NOTCH and WNT expression 
colData(salmon_quants_hg38.noimmune)$notch <- assays(salmon_quants_hg38.noimmune)$logcounts[rowData(salmon_quants_hg38.noimmune)$symbol == "NOTCH1",]
colData(salmon_quants_hg38.hiseq.noimmune)$notch <- assays(salmon_quants_hg38.hiseq.noimmune)$logcounts[rowData(salmon_quants_hg38.hiseq.noimmune)$symbol == "NOTCH1",]

colData(salmon_quants_hg38.noimmune)$axin2 <- assays(salmon_quants_hg38.noimmune)$logcounts[rowData(salmon_quants_hg38.noimmune)$symbol == "AXIN2",]
colData(salmon_quants_hg38.hiseq.noimmune)$axin2 <- assays(salmon_quants_hg38.hiseq.noimmune)$logcounts[rowData(salmon_quants_hg38.hiseq.noimmune)$symbol == "AXIN2",]

colData(salmon_quants_hg38.noimmune)$hes1 <- assays(salmon_quants_hg38.noimmune)$logcounts[rowData(salmon_quants_hg38.noimmune)$symbol == "HES1",]
colData(salmon_quants_hg38.hiseq.noimmune)$hes1 <- assays(salmon_quants_hg38.hiseq.noimmune)$logcounts[rowData(salmon_quants_hg38.hiseq.noimmune)$symbol == "HES1",]

to_plot.nova <- data.frame(pc1 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,1],
                           pc2 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,2],
                           pseudotime = sling_time.nova,
                           clusters = salmon_quants_hg38.noimmune$cluster.collapse)

to_plot.hiseq <- data.frame(pc1 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,1],
                            pc2 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,2],
                            pseudotime = sling_time.hiseq,
                            clusters = salmon_quants_hg38.hiseq.noimmune$cluster.collapse)

# plot
library(ggplot2)
library(viridis)

# nova-seq pat 1
ggplot(to_plot.nova, aes(x = pc1, y = pc2, color = pseudotime)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "inferno",
                        limits = c(0,1),
                        breaks = seq(0,1,1)) +
  geom_path(
    data = curve_df.nova,
    aes(x = Dim1, y = Dim2, group = lineage),
    color = "black",
    linewidth = 1.5
  ) +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  )

# hiseq pat 2
ggplot(to_plot.hiseq, aes(x = pc1, y = pc2, color = pseudotime)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "inferno",
                        limits = c(0,1),
                        breaks = seq(0,1,1)) +
  geom_path(
    data = curve_df.hiseq,
    aes(x = Dim1, y = Dim2, group = lineage),
    color = "black",
    linewidth = 1.5
  ) +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  )

to_plot.nova <- data.frame(pc1 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,1],
                           pc2 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,2],
                           pseudotime = sling_time.nova,
                           clusters = salmon_quants_hg38.noimmune$cluster)

to_plot.hiseq <- data.frame(pc1 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,1],
                            pc2 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,2],
                            pseudotime = sling_time.hiseq,
                            clusters = salmon_quants_hg38.hiseq.noimmune$cluster)

save(to_plot.nova, to_plot.hiseq,
     file = "~/Documents/manuscripts/storm_seq/fte_analysis/principle_curves_pat1_pat2_with_pseudotime_20250616.rda")

ggplot(to_plot.nova, aes(x = pc1, y = pc2, color = clusters)) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "None",
    panel.grid = element_blank()
  )

ggplot(to_plot.hiseq, aes(x = pc1, y = pc2, color = clusters)) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "None",
    panel.grid = element_blank()
  )

# PCA plots of NOTCH/WNT signaling
to_plot.nova <- data.frame(pc1 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,1],
                           pc2 = reducedDim(salmon_quants_hg38.noimmune, "PCAsub")[,2],
                           notch = salmon_quants_hg38.noimmune$notch,
                           axin2 = salmon_quants_hg38.noimmune$axin2)

to_plot.hiseq <- data.frame(pc1 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,1],
                            pc2 = reducedDim(salmon_quants_hg38.hiseq.noimmune, "PCAsub")[,2],
                            notch = salmon_quants_hg38.hiseq.noimmune$notch,
                            axin2 = salmon_quants_hg38.hiseq.noimmune$axin2)
# Notch1
ggplot(to_plot.nova, aes(x = pc1, y = pc2, color = notch)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "log2(NOTCH1)") +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

ggplot(to_plot.hiseq, aes(x = pc1, y = pc2, color = notch)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "log2(NOTCH1)") +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

# Axin2
ggplot(to_plot.nova, aes(x = pc1, y = pc2, color = axin2)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "log2(AXIN2)",
                        breaks = c(0.0, 4.0, 8.0, 12.0),
                        limits = c(0.0, 12.0),
                        labels = scales::number_format(accuracy = 0.1)) +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

ggplot(to_plot.hiseq, aes(x = pc1, y = pc2, color = axin2)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "log2(AXIN2)",
                        breaks = c(0.0, 4.0, 8.0, 12.0),
                        limits = c(0.0, 12.0),
                        labels = scales::number_format(accuracy = 0.1)) +
  #ggtitle("STORM-seq FTE PCA") +
  xlab("PC 1") +
  ylab("PC 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

# let's just go straight to gsea
# Packages

library(matrixStats)
library(msigdbr)
library(GSVA)
library(BiocParallel)
library(fgsea)
library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

## ---- inputs ----
sce <- salmon_quants_hg38
cluster_col <- "cluster"  # <-- change if your cluster column has a different name

stopifnot("logcounts" %in% assayNames(sce))
stopifnot("symbol" %in% colnames(rowData(sce)))
stopifnot(cluster_col %in% colnames(colData(sce)))

## ---- collapse Ensembl -> unique SYMBOL rows (per-cell column-wise max) ----
mat <- as.matrix(assay(sce, "logcounts"))    # genes x cells (log-scale)
sym <- rowData(sce)$symbol
keep <- !is.na(sym) & sym != ""
mat  <- mat[keep, , drop = FALSE]
sym  <- sym[keep]

collapse_by_symbol <- function(M, syms){
  idx_list <- split(seq_along(syms), syms)   # list of row indices per SYMBOL
  # For each symbol's rows, take column-wise max (fast):
  out <- vapply(idx_list, function(ix) matrixStats::colMaxs(M[ix, , drop = FALSE]),
                numeric(ncol(M)))
  out <- t(out)                               # SYMBOL x cells
  rownames(out) <- names(idx_list)
  colnames(out) <- colnames(M)
  out
}
mat_sym <- collapse_by_symbol(mat, sym)       # SYMBOL x cells

## ---- cluster pseudobulk mean (genes x clusters) ----
clusters <- factor(colData(sce)[[cluster_col]])  # keep level order stable
split_idx <- split(seq_along(clusters), clusters, drop = FALSE)
avg_by_cluster <- vapply(split_idx, function(ix) Matrix::rowMeans(mat_sym[, ix, drop = FALSE]),
                         numeric(nrow(mat_sym)))
avg_by_cluster <- as.matrix(avg_by_cluster)      # genes x clusters
rownames(avg_by_cluster) <- rownames(mat_sym)

## ---- gene sets from msigdbr (restrict to your picks) ----
wanted_stem <- c(
  "WONG_EMBRYONIC_STEM_CELL_CORE",
  "BENPORATH_ES_1",
  "BENPORATH_ES_CORE_NINE",
  "MALTA_CURATED_STEMNESS_MARKERS",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS",
  "WONG_ADULT_TISSUE_STEM_MODULE",
  "GAO_LARGE_INTESTINE_ADULT_CE_OLFM4HIGH_STEM_CELL"
)

# Ciliated cell programs (cilia/axoneme/beat)
wanted_cilia <- c(
  "GOBP_CILIUM_MOVEMENT",
  "GOBP_CILIUM_ORGANIZATION",
  "GOBP_CILIUM_OR_FLAGELLUM_DEPENDENT_CELL_MOTILITY",
  "CILIARY_BASAL_BODY_ORGANIZATION"
)

# Secretory programs (secretion/vesicle trafficking)
wanted_secretory <- c(
  "GOBP_REGULATED_EXOCYTOSIS",
  "GOBP_SNARE_COMPLEX_ASSEMBLY",
  "GOBP_APICAL_PLASMA_MEMBRANE_PROTEIN_LOCALIZATION",
  "GOBP_PROTEIN_O_LINKED_GLYCOSYLATION",
  "GOBP_PROTEIN_N_LINKED_GLYCOSYLATION",
  "GOBP_MUCIN_BIOSYNTHETIC_PROCESS",
  "GOBP_SECRETORY_GRANULE_ORGANIZATION"
)

wanted_all <- c(wanted_stem, wanted_cilia, wanted_secretory)

msig_all <- msigdbr(species = "Homo sapiens") |>
  dplyr::filter(gs_name %in% wanted_all) |>
  dplyr::distinct(gs_name, gene_symbol)

gsets <- split(msig_all$gene_symbol, msig_all$gs_name)

## ---- optional: add your FT mini-set (from the FT organoid paper) ----
gsets[["KESSLER_FT_STEMNESS_ORGANOID"]] <- c("OLFM4","LGR6","HES1","AXIN2","TERT","RNF43","FZD2","FZD7","EPHA4","SMO")
gsets[["FT_SECRETORY_MINI"]] <- c("OVGP1","PAX8","SLPI","LCN2","WFDC2","GPX3","MUC1","MUC16","KRT7","MSLN")
## ---- intersect gene sets with gene universe & enforce min size ----
gene_universe <- rownames(avg_by_cluster)
gsets_eff <- lapply(gsets, intersect, gene_universe)
min_size <- 10
gsets_eff <- gsets_eff[vapply(gsets_eff, length, integer(1)) >= min_size]

## ---- (optional) prefilter constant genes to speed GSVA and reduce warnings ----
keep_var <- matrixStats::rowSds(avg_by_cluster) > 0
avg_by_cluster_f <- avg_by_cluster[keep_var, , drop = FALSE]
# Recompute gene sets against filtered universe
gene_universe_f <- rownames(avg_by_cluster_f)
gsets_eff <- lapply(gsets_eff, intersect, gene_universe_f)
gsets_eff <- gsets_eff[vapply(gsets_eff, length, integer(1)) >= min_size]

## ---- ssGSEA with new API (ssgseaParam -> gsva) ----
param_ssgsea <- ssgseaParam(
  exprData  = avg_by_cluster_f,  # genes x clusters
  geneSets  = gsets_eff,
  minSize   = min_size,
  maxSize   = 5000,
  alpha     = 0.25,              # ssGSEA exponent (Barbie et al.)
  normalize = TRUE               # normalize by (max - min)
)

ssgsea_scores <- gsva(
  param   = param_ssgsea,
  verbose = TRUE,
  BPPARAM = SerialParam(progressbar = TRUE)   # use MulticoreParam(workers=8) if desired
)
# ssgsea_scores: genesets x clusters (matrix); attr(,"geneSets") holds the sets

## ---- tidy & plot ssGSEA heatmap ----
df_ss <- as.data.frame(t(ssgsea_scores)) |>
  tibble::rownames_to_column("cluster") |>
  tidyr::pivot_longer(-cluster, names_to = "geneset", values_to = "ssgsea")

df_ss_z <- df_ss %>%
  group_by(geneset) %>%
  mutate(zscore = (ssgsea - mean(ssgsea)) / sd(ssgsea)) %>%
  ungroup()

ggplot(df_ss, aes(cluster, geneset, fill = ssgsea)) +
  geom_tile() +
  scale_fill_viridis_c(name = "ssGSEA") +
  theme_minimal(base_size = 12) +
  labs(title = "ssGSEA (cluster pseudobulk)")

# plot heatmap with z-scores
ggplot(df_ss_z, aes(cluster, geneset, fill = zscore)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "magma",
    name   = "ssGSEA (z-score)",
    limits = c(-2, 2), oob = scales::squish) +
  labs(title = "ssGSEA") +
  xlab("Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
    )


## ---- inputs ----
sce <- salmon_quants_hg38.hiseq
cluster_col <- "cluster"  # <-- change if your cluster column has a different name

stopifnot("logcounts" %in% assayNames(sce))
stopifnot("symbol" %in% colnames(rowData(sce)))
stopifnot(cluster_col %in% colnames(colData(sce)))

## ---- collapse Ensembl -> unique SYMBOL rows (per-cell column-wise max) ----
mat <- as.matrix(assay(sce, "logcounts"))    # genes x cells (log-scale)
sym <- rowData(sce)$symbol
keep <- !is.na(sym) & sym != ""
mat  <- mat[keep, , drop = FALSE]
sym  <- sym[keep]

collapse_by_symbol <- function(M, syms){
  idx_list <- split(seq_along(syms), syms)   # list of row indices per SYMBOL
  # For each symbol's rows, take column-wise max (fast):
  out <- vapply(idx_list, function(ix) matrixStats::colMaxs(M[ix, , drop = FALSE]),
                numeric(ncol(M)))
  out <- t(out)                               # SYMBOL x cells
  rownames(out) <- names(idx_list)
  colnames(out) <- colnames(M)
  out
}
mat_sym <- collapse_by_symbol(mat, sym)       # SYMBOL x cells

## ---- cluster pseudobulk mean (genes x clusters) ----
clusters <- factor(colData(sce)[[cluster_col]])  # keep level order stable
split_idx <- split(seq_along(clusters), clusters, drop = FALSE)
avg_by_cluster <- vapply(split_idx, function(ix) Matrix::rowMeans(mat_sym[, ix, drop = FALSE]),
                         numeric(nrow(mat_sym)))
avg_by_cluster <- as.matrix(avg_by_cluster)      # genes x clusters
rownames(avg_by_cluster) <- rownames(mat_sym)

## ---- gene sets from msigdbr (restrict to your picks) ----
wanted <- c(
  "WONG_EMBRYONIC_STEM_CELL_CORE",
  "BENPORATH_ES_1",
  "BENPORATH_ES_CORE_NINE",
  "MALTA_CURATED_STEMNESS_MARKERS",
  "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS",
  "WONG_ADULT_TISSUE_STEM_MODULE",
  "GAO_LARGE_INTESTINE_ADULT_CE_OLFM4HIGH_STEM_CELL"
)

msig_all <- msigdbr(species = "Homo sapiens") |>
  dplyr::filter(gs_name %in% wanted) |>
  dplyr::distinct(gs_name, gene_symbol)

gsets <- split(msig_all$gene_symbol, msig_all$gs_name)

## ---- optional: add your FT mini-set (from the FT organoid paper) ----
gsets[["KESSLER_FT_STEMNESS_ORGANOID"]] <- c("OLFM4","LGR6","HES1","AXIN2","TERT","RNF43","FZD2","FZD7","EPHA4","SMO")

## ---- intersect gene sets with gene universe & enforce min size ----
gene_universe <- rownames(avg_by_cluster)
gsets_eff <- lapply(gsets, intersect, gene_universe)
min_size <- 10
gsets_eff <- gsets_eff[vapply(gsets_eff, length, integer(1)) >= min_size]

## ---- (optional) prefilter constant genes to speed GSVA and reduce warnings ----
keep_var <- matrixStats::rowSds(avg_by_cluster) > 0
avg_by_cluster_f <- avg_by_cluster[keep_var, , drop = FALSE]
# Recompute gene sets against filtered universe
gene_universe_f <- rownames(avg_by_cluster_f)
gsets_eff <- lapply(gsets_eff, intersect, gene_universe_f)
gsets_eff <- gsets_eff[vapply(gsets_eff, length, integer(1)) >= min_size]

## ---- ssGSEA with new API (ssgseaParam -> gsva) ----
param_ssgsea <- ssgseaParam(
  exprData  = avg_by_cluster_f,  # genes x clusters
  geneSets  = gsets_eff,
  minSize   = min_size,
  maxSize   = 5000,
  alpha     = 0.25,              # ssGSEA exponent (Barbie et al.)
  normalize = TRUE               # normalize by (max - min)
)

ssgsea_scores <- gsva(
  param   = param_ssgsea,
  verbose = TRUE,
  BPPARAM = SerialParam(progressbar = TRUE)   # use MulticoreParam(workers=8) if desired
)
# ssgsea_scores: genesets x clusters (matrix); attr(,"geneSets") holds the sets

## ---- tidy & plot ssGSEA heatmap ----
df_ss <- as.data.frame(t(ssgsea_scores)) |>
  tibble::rownames_to_column("cluster") |>
  tidyr::pivot_longer(-cluster, names_to = "geneset", values_to = "ssgsea")

df_ss_z <- df_ss %>%
  group_by(geneset) %>%
  mutate(zscore = (ssgsea - mean(ssgsea)) / sd(ssgsea)) %>%
  ungroup()

ggplot(df_ss, aes(cluster, geneset, fill = ssgsea)) +
  geom_tile() +
  scale_fill_viridis_c(name = "ssGSEA") +
  theme_minimal(base_size = 12) +
  labs(title = "ssGSEA (cluster pseudobulk)")

# plot heatmap with z-scores
ggplot(df_ss_z, aes(cluster, geneset, fill = zscore)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "magma",
    name   = "ssGSEA (z-score)",
    limits = c(-2, 2), oob = scales::squish) +
  labs(title = "ssGSEA") +
  xlab("Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
