# ingest and process the FTE TE results

# patient 1
# mapped using Ensembl and repeatmasker with Ayush's TE pipeline
te_cpm_filt <- read.delim("count_matrix_filtered_cpm_TE.tsv.gz",
                          header = TRUE)
rownames(te_cpm_filt) <- te_cpm_filt$Geneid

# quick filter all zero entries
# there are 306 cells after QC and filtering...
table(rowSums(te_cpm_filt[,c(1:306)]) > 0)
# FALSE    TRUE 
# 1789937 2748141 
te_cpm_filt <- te_cpm_filt[rowSums(te_cpm_filt[,c(1:306)]) > 0,]

# turn this into a SingleCellExperiment
library(SingleCellExperiment)
library(Matrix)
te_cpm_filt.sce <- SingleCellExperiment(assay = list(counts = te_cpm_filt[,c(1:306)]))
assay(te_cpm_filt.sce, "logcounts") <- log2(counts(te_cpm_filt.sce) + 1) 
my_rowdata <- data.frame(Geneid = te_cpm_filt$Geneid,
                         repName = te_cpm_filt$repName,
                         repClass = te_cpm_filt$repClass,
                         repFamily = te_cpm_filt$repFamily)
rowData(te_cpm_filt.sce) <- my_rowdata
saveRDS(te_cpm_filt.sce, file = "pat1_novaseq_fte_te_cpm_sce_ens101_raw.rds")
# te_cpm_filt.sce <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/pat1_novaseq_fte_te_cpm_sce_ens101_raw.rds")
# filter again on HPC since this is a large object, even after filtering for non-zero entries
# subset to just LINEs, LTRs, and SINES
# te_cpm_filt.sce <- te_cpm_filt.sce[rowData(te_cpm_filt.sce)$repClass %in% c("LINE",
#                                                                             "LTR",
#                                                                             "SINE"),]

# filter TEs to those that are expressed in at least 10% of cells
# keep_genes <- rowMeans(logcounts(te_cpm_filt.sce) > 0) >= 0.1
# FALSE    TRUE
# 2402144   67483

# filter
# te_cpm_filt.sce <- te_cpm_filt.sce[keep_genes, ]

te_cpm_filt.sce <- readRDS("pat1_novaseq_fte_te_cpm_sce_ens101_ten_perc_filt_line_sine_ltr.rds")

# now we need to drag in the filtered patient 1 expression object
# so that we can add clusters and filtered cells
pat1_sce <- readRDS("cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")
colnames(pat1_sce) <- gsub("_quant", "", colnames(pat1_sce))

# filter to cells found in the TE sce
pat1_sce.filt_vec <- colnames(pat1_sce) %in% colnames(te_cpm_filt.sce)
pat1_sce <- pat1_sce[,pat1_sce.filt_vec]

# match and filter cells
pat1_sce.m <- match(colnames(te_cpm_filt.sce),
                    colnames(pat1_sce))
pat1_sce <- pat1_sce[,pat1_sce.m]

# careful
all(colnames(te_cpm_filt.sce) == colnames(pat1_sce))
# TRUE

# tack on the cluster IDs from the pat1_sce
te_cpm_filt.sce$cluster <- pat1_sce$cluster

# look for marker TEs
library(scran)
te_markers <- findMarkers(te_cpm_filt.sce,
                          groups = te_cpm_filt.sce$cluster,
                          test.type = "wilcox",
                          pval.type = "some",
                          assay.type = "logcounts")

# extract significant TE markers
te_markers.sig <- unique(unlist(lapply(te_markers, function(x) {
  x.sig <- x[x$FDR < 0.05 & x$summary.AUC >= 0.9,]
  if (nrow(x.sig) > 0) {
    return(rownames(x.sig))
  } else {
    return(NA)
  }
})))

# find a strong example for ciliated cells?
# this can be "relaxed" by removing the AUC filter above and
# filter by AUC >= 0.8 here
te_markers.cil <- te_markers[names(te_markers) %in% c("1", "4")]
te_markers.cil <- lapply(te_markers.cil, function(x) {
  return(x[x$FDR < 0.05 & x$summary.AUC >= 0.9,])
})

# find a strong example for secretory cells?
te_markers.sec <- te_markers[names(te_markers) %in% c("2", "7", "8")]
te_markers.sec <- lapply(te_markers.sec, function(x) {
  return(x[x$FDR < 0.05 & x$summary.AUC >= 0.8,])
})

# filter to intergenic TEs here too
# note: this is based on Ensembl 101 annotations with repeat masker
te_annot <- read.delim("intergenic_intronic_tes.txt")
te_annot.sig_te <- te_annot[te_annot$name %in% te_markers.sig,]
te_annot.sig_te.intergenic <- te_annot.sig_te[te_annot.sig_te$annot %in% "intergenic",]

te_markers.cil.intergenic <- lapply(te_markers.cil, function(x) {
  return(x[rownames(x) %in% te_annot.sig_te.intergenic$name,])
})

te_markers.sig.intergenic <- lapply(te_markers, function(x) {
  return(x[rownames(x) %in% te_annot.sig_te.intergenic$name,])
})

# make a heatmap
library(scater)
# ciliated
scater::plotExpression(te_cpm_filt.sce, features = "te_3881379",
                       colour_by = "cluster", x = "cluster")
scater::plotExpression(te_cpm_filt.sce, features = "te_1075180",
                       colour_by = "cluster", x = "cluster")

# secretory
scater::plotExpression(te_cpm_filt.sce, features = "te_2160307",
                       colour_by = "cluster", x = "cluster")
scater::plotExpression(te_cpm_filt.sce, features = "te_2155624",
                       colour_by = "cluster", x = "cluster")


reducedDim(te_cpm_filt.sce, "densMAP") <- reducedDim(pat1_sce, "densMAP")

library(ggplot2)
dmap <- reducedDim(te_cpm_filt.sce, "densMAP")
to_plot <- data.frame(dim1 = dmap[,1],
                      dim2 = dmap[,2])

# gather the ciliated cell marker TE expression
te_cil <- unique(unlist(lapply(te_markers.cil.intergenic, function(x) rownames(x))))
# te_annot.cil_markers.relax <- te_annot.sig_te.intergenic[te_annot.sig_te.intergenic$name %in% te_cil,]

te_cil.logcounts <- logcounts(te_cpm_filt.sce)[rownames(logcounts(te_cpm_filt.sce)) %in% te_cil,]
to_plot <- cbind(to_plot,
                 t(te_cil.logcounts))

# remove the immune pop for plotting
immune <- colnames(te_cpm_filt.sce)[te_cpm_filt.sce$cluster %in% "6"]
to_plot <- to_plot[!rownames(to_plot) %in% immune,]

# te_3881379
# te_1075180
ggplot(to_plot, aes(x = dim1, y = dim2, color = te_1075180)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  ggtitle("L1PA3 (LINE)") +
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    #legend.position = "None",
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

# heatmap of TE specific expression
library(ComplexHeatmap)

te_sig.intergenic <- unique(unlist(lapply(te_markers.sig.intergenic, function(x) rownames(x))))
sig_fte_te_exp <- logcounts(te_cpm_filt.sce)[rownames(te_cpm_filt.sce) %in% te_sig.intergenic,]

getJetColors <- function(circlize.cols=NULL) {
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (!is.null(circlize.cols)) {
    jet.colors <- circlize::colorRamp2(seq(0,10,length.out = 9), c("#00007F", "blue", "#007FFF", "cyan",
                                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  return(jet.colors)
}

# remove the immune cluster
sig_fte_te_exp <- sig_fte_te_exp[,!te_cpm_filt.sce$cluster %in% "6",]
noimmune_clust <- te_cpm_filt.sce$cluster[!te_cpm_filt.sce$cluster %in% "6"]

heatmap_legend_params <- list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 14),
                              legend_width = unit(5, "cm"),
                              legend_height = unit(5, "cm"))

Heatmap(as.matrix(sig_fte_te_exp),
        column_split = noimmune_clust,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params) 

# write out the per-cluster marker TEs from patient 1
pat1_path <- "~/my/path"
lapply(names(te_markers.sig.intergenic), function(x) {
  marker.sub <- te_markers.sig.intergenic[[x]]
  write.table(marker.sub,
              file = paste0(pat1_path, "pat1_fte_te_markers_cluster_", x, ".txt"),
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE,
              sep = '\t')
  message("Written: ",
          paste0(pat1_path, "pat1_fte_te_markers_cluster_", x, ".txt"))
})

# write out the unique TE "markers"
write.table(te_sig.intergenic,
            file = "pat1_sig_te_markers_across_clusters.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

# look in the second patient too
# patient 2
# mapped using Ensembl and repeatmasker with Ayush's TE pipeline
te_cpm_filt2 <- read.delim("count_matrix_filtered_cpm_TE.tsv",
                           header = TRUE)
rownames(te_cpm_filt2) <- te_cpm_filt2$Geneid

# quick filter all zero entries
# again, 335 cells after QC and filtering...
table(rowSums(te_cpm_filt2[,c(1:335)]) > 0)
# FALSE    TRUE
# 2621692 1916386
te_cpm_filt2 <- te_cpm_filt2[rowSums(te_cpm_filt2[,c(1:335)]) > 0,]

# turn this into a SingleCellExperiment
library(SingleCellExperiment)
library(Matrix)
te_cpm_filt2.sce <- SingleCellExperiment(assay = list(counts = te_cpm_filt2[,c(1:335)]))
assay(te_cpm_filt2.sce, "logcounts") <- log2(counts(te_cpm_filt2.sce) + 1) 
my_rowdata2 <- data.frame(Geneid = te_cpm_filt2$Geneid,
                         repName = te_cpm_filt2$repName,
                         repClass = te_cpm_filt2$repClass,
                         repFamily = te_cpm_filt2$repFamily)
rowData(te_cpm_filt2.sce) <- my_rowdata2
saveRDS(te_cpm_filt2.sce, file = "pat2_hiseq_fte_te_cpm_sce_ens101_raw.rds")

# drag in the sce
te_cpm_filt2.sce <- readRDS("pat2_hiseq_fte_te_cpm_sce_ens101_raw.rds")

# filter to known markers
te_cpm_filt2.sce.sig <- te_cpm_filt2.sce[rownames(te_cpm_filt2.sce) %in% te_sig.intergenic,]
# grabbed all but _one_ marker
# class: SingleCellExperiment
# dim: 248 335
# metadata(0):
#   assays(2): counts logcounts
# rownames(248): te_36433 te_76665 ... te_4636987 te_4758258
# rowData names(4): Geneid repName repClass repFamily
# colnames(335): SE6052_SA56912_S1_L001 SE6052_SA56914_S3_L001 ...
# SE6054_SA56910_S191_L002 SE6054_SA56911_S192_L002
# colData names(0):
#   reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

# pull in the pat 2 gene expression sce to add the cluster level info
pat2_sce <- readRDS("cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")
colnames(pat2_sce) <- gsub(".merged_quant|_quant", "", colnames(pat2_sce))
colnames(pat2_sce) <- gsub("_R1_001", "", colnames(pat2_sce))

# filter first since there are fewer cells (QC'd already)
# in the pat2_sce than in the te sce
te_cpm_filt2.sce.sig <- te_cpm_filt2.sce.sig[,colnames(te_cpm_filt2.sce.sig) %in% colnames(pat2_sce)]

pat2_sce.m <- match(colnames(te_cpm_filt2.sce.sig),
                    colnames(pat2_sce))
pat2_sce <- pat2_sce[,pat2_sce.m]

# careful
all(colnames(pat2_sce) == colnames(te_cpm_filt2.sce.sig))
# TRUE

# tack on the clusters and embeddings
te_cpm_filt2.sce.sig$cluster <- pat2_sce$cluster
reducedDim(te_cpm_filt2.sce.sig, "TSNE") <- reducedDim(pat2_sce, "TSNE")
reducedDim(te_cpm_filt2.sce.sig, "denSNE") <- reducedDim(pat2_sce, "denSNE")

# save
saveRDS(te_cpm_filt2.sce.sig, "pat2_hiseq_fte_te_cpm_sce_ens101_te_markers_with_embeddings.rds")

library(ComplexHeatmap)

getJetColors <- function(circlize.cols=NULL) {
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (!is.null(circlize.cols)) {
    jet.colors <- circlize::colorRamp2(seq(0,10,length.out = 9), c("#00007F", "blue", "#007FFF", "cyan",
                                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  return(jet.colors)
}

heatmap_legend_params <- list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 14),
                              legend_width = unit(5, "cm"),
                              legend_height = unit(5, "cm"))

Heatmap(as.matrix(logcounts(te_cpm_filt2.sce.sig)),
        column_split = te_cpm_filt2.sce.sig$cluster,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params)

library(ggplot2)
dmap <- reducedDim(te_cpm_filt2.sce.sig, "TSNE")
to_plot <- data.frame(dim1 = dmap[,1],
                      dim2 = dmap[,2])

te_pat2.logcounts <- logcounts(te_cpm_filt2.sce.sig)
to_plot <- cbind(to_plot,
                 t(te_pat2.logcounts))

# te_3881379
# te_1075180
ggplot(to_plot, aes(x = dim1, y = dim2, color = te_3881379)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  ggtitle("L1PA3 (LINE)") +
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    #legend.position = "None",
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

# fed the relaxed ciliated TE markers to FIMO and parse the results
fimo <- read.delim("fimo.tsv")

# filter to significant hits
fimo.sig <- fimo[fimo$q.value < 0.05,]

# looks like zinc-finger proteins dominate 
