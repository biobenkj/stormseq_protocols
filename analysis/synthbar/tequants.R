## process the K562 TE data for synthbar manuscript

## modify the config file for using the resulting STARsolo bams
rnames <- read.delim("read_names.tsv", header = TRUE)
# [1] "cells"           "paired"          "raw_reads_1"     "raw_reads_2"
# [5] "trim"            "trimmed_reads_1" "trimmed_reads_2" "align"
# [9] "aligned_files"   "fc"              "fc_files"        "sc"

# gather available cells
cell_bams <- list.files(path = "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/k562",
                        pattern = "*.out.bam$")
cells <- gsub("_soloAligned.sortedByCoord.out.bam",
              "",
              cell_bams)
rnames_new <- data.frame(cells = cells,
                         paired = rep("TRUE",
                                      length(cells)),
                         raw_reads_1 = rep("NULL",
                                           length(cells)),
                         raw_reads_2 = rep("NULL",
                                           length(cells)),
                         trim = rep("FALSE",
                                    length(cells)),
                         trimmed_reads_1 = rep("NULL",
                                               length(cells)),
                         trimmed_reads_2 = rep("NULL",
                                               length(cells)),
                         align = rep("FALSE",
                                     length(cells)),
                         aligned_files = cell_bams,
                         fc = rep("TRUE",
                                  length(cells)),
                         fc_files = rep("NULL",
                                        length(cells)),
                         sc = rep("TRUE",
                                  length(cells)))
write.table(rnames_new,
            file = "read_names.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

## run the pipeline
source("TE_quantification_pipeline/TE_quantification_pipeline.R")
target_dir <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/k562/storm_k562_tequant"
config_file <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/k562/TE_quantification_pipeline/supporting_files/config_file.tsv"
read_names_file <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/k562/TE_quantification_pipeline/supporting_files/read_names.tsv"
sample <- "storm_k562"
create_count_matrix <- "TRUE"
scuttle_filter <- "TRUE"
log_enrich <- "TRUE"
quantify_TEs(config_file, read_names_file, sample, target_dir, create_count_matrix, scuttle_filter, log_enrich)

## import filtered CPM
k562_te_quants <- read.delim("~/Documents/manuscripts/synthbar_2025/k562_tequants/count_matrix_filtered_cpm_TE.tsv.gz",
                             header = TRUE)

## adding genes as rownames
rownames(k562_te_quants) <- k562_te_quants$Geneid

# quick filter all zero entries
table(rowSums(k562_te_quants[,c(1:95)]) > 0)
# FALSE    TRUE 
# 4078846  459232 
k562_te_quants <- k562_te_quants[rowSums(k562_te_quants[,c(1:95)]) > 0,]

# filter to TEs in at least 20% of cells
keep_genes <- rowMeans(k562_te_quants[,c(1:95)] > 0) >= 0.2
# FALSE   TRUE 
# 453430   5802

k562_te_quants.filt <- k562_te_quants[keep_genes,]

# subset to LINE, LTR, and SINEs
# distributions
# DNA       LINE        LTR         RC Retroposon       SINE 
# 381       2540       1588          2         53       1238

k562_te_quants.filt.line <- k562_te_quants.filt[k562_te_quants.filt$repClass %in% "LINE",]
k562_te_quants.filt.ltr <- k562_te_quants.filt[k562_te_quants.filt$repClass %in% "LTR",]
k562_te_quants.filt.sine <- k562_te_quants.filt[k562_te_quants.filt$repClass %in% "SINE",]

# log it
k562_te_quants.filt.line.log <- log2(k562_te_quants.filt.line[,1:95] + 1)
k562_te_quants.filt.ltr.log <- log2(k562_te_quants.filt.ltr[,1:95] + 1)
k562_te_quants.filt.sine.log <- log2(k562_te_quants.filt.sine[,1:95] + 1)

# combine and write out table
k562_tes_out <- rbind(k562_te_quants.filt.line.log,
                      k562_te_quants.filt.ltr.log,
                      k562_te_quants.filt.sine.log)

## plot as a heatmap
library(ComplexHeatmap)

getJetColors <- function(circlize.cols=NULL) {
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (!is.null(circlize.cols)) {
    jet.colors <- circlize::colorRamp2(seq(0,6,length.out = 9), c("#00007F", "blue", "#007FFF", "cyan",
                                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  return(jet.colors)
}

heatmap_legend_params <- list(title_gp = gpar(fontsize = 14, fontface = "bold"),
                              labels_gp = gpar(fontsize = 14),
                              legend_width = unit(3, "cm"),
                              legend_height = unit(3, "cm"))

# LINE
Heatmap(as.matrix(k562_te_quants.filt.line.log),
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params)

# LTR
Heatmap(as.matrix(k562_te_quants.filt.ltr.log),
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params)

# SINE
Heatmap(as.matrix(k562_te_quants.filt.sine.log),
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params)



