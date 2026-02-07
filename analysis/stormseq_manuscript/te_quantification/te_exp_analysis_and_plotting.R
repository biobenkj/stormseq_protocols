## process TE data

## TEs UMI STORM from STARsolo
library(SingleCellExperiment)
library(ggplot2)
library(scales)
library(data.table)
library(stringr)

# filtered_log_enrichment_STORM <- fread("/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/storm_te_quant/filtered_log_enrichment.tsv", header = TRUE)
# 
# filtered_log_enrichment_STORM$repClass.x <- factor(filtered_log_enrichment_STORM$repClass.x, 
#                                                    levels = str_sort(unique(filtered_log_enrichment_STORM$repClass.x)))
# 
# filtered_log_enrichment_STORM_plot <- ggplot(filtered_log_enrichment_STORM, aes(repFamily, cells)) + 
#   geom_tile(aes(fill = log_enrichment)) +
#   facet_grid(~ repClass.x, 
#              scales = "free", space = "free") +
#   scale_fill_gradientn(name = "log2(Enrichment)", colours = c("dodgerblue4","whitesmoke","firebrick3"), 
#                        values = rescale(c(-3,0,3)),
#                        guide = "colorbar", limits=c(-3,3)) +
#   ggtitle("STORM-seq") +
#   theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
#         strip.background = element_blank(),
#         strip.text.x = element_text(size=12, face="bold"),
#         panel.border = element_rect(color = "black", fill = NA, size = 1), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), 
#         axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
#         axis.text.y=element_blank()) +
#   xlab("Repeat Family") +
#   ylab("Condition")
# filtered_log_enrichment_STORM_plot

## code snippet from Josh
library(ggplot2) 
library(RColorBrewer) 
library(plyr)
library(ComplexHeatmap) 
library(circlize) 
library(reshape2) 
library(DESeq2)
library(matrixStats) 

mypalette <- brewer.pal(12,"Paired")
mypalette2 <- brewer.pal(12,"Set3")
mypalette3 <- brewer.pal(8,"Pastel1")
mypalette4 <- brewer.pal(8,"Set2")
mypalette5 <- brewer.pal(11,"Spectral")
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }


## need to import dataSTORM and filter as appropriate
storm_sce <- readRDS("~/Documents/kb_python_storm_test/sce_filt_ens101_noround.rds")
#storm_sce <- storm_sce[,storm_sce$CellType %in% "K562"]
storm_sce <- storm_sce[,storm_sce$CellType %in% "HEK293T"]
## re-filter to only cells that passed QC threshold for gene expression
dataSTORM.filt <- dataSTORM[,colnames(dataSTORM)[1:343] %in% colnames(storm_sce)]
## some reason we lost a couple columns
# this is fine for some reason?
#dataSTORM.filt[,105:106] <- NULL
dataSTORM.filt <- cbind(dataSTORM.filt,
                        dataSTORM[,344:347])

#BulkTotal <- read.table("/Volumes/projects/shen/projects/2022_08_01_Josh_Ben_Ayush_STORM_TE_project/K562_data/bulk/count_matrix_cpm_TE.tsv",sep ="\t", header = T, stringsAsFactors = F)
BulkTotal <- read.table("/Volumes/projects/shen/projects/2022_08_01_Josh_Ben_Ayush_STORM_TE_project/HEK293_data/HEK293_total_bulk/count_matrix_cpm_TE.tsv",sep ="\t", header = T, stringsAsFactors = F)


## other technologies
# ss2 <- read.table("/Volumes/projects/shen/projects/2022_08_01_Josh_Ben_Ayush_STORM_TE_project/K562_data/Natarajan/count_matrix_filtered_cpm_TE.tsv", sep = '\t', header = T, stringsAsFactors = F)
# vasa <- read.table("~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/vasa/te_quants/count_matrix_filtered_cpm_TE.tsv", sep = '\t', header = T, stringsAsFactors = F)
sstotal <- read.table("/Users/benjamin.johnson/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/sstotal/te_quants/count_matrix_filtered_cpm_TE.tsv", sep = '\t', header = T, stringsAsFactors = F)


## subset vasa to just the hek cells
# annotate cell types and controls
# note: P1 and P2 have same controls

colnames(vasa)[1:676] <- c(paste0(colnames(vasa)[1:165], "_K562"),
                           paste0(colnames(vasa)[166:339], "_HEK293T"),
                           paste0(colnames(vasa)[340:499], "_K562"),
                           paste0(colnames(vasa)[500:676], "_HEK293T"))

p1_neg_controls <- c("VAI.MR.v001_HNVCGBGXN_S3_009_K562", "VAI.MR.v001_HNVCGBGXN_S3_047_K562",
                     "VAI.MR.v001_HNVCGBGXN_S3_108_K562", "VAI.MR.v001_HNVCGBGXN_S3_123_K562",
                     "VAI.MR.v001_HNVCGBGXN_S3_207_HEK293T", "VAI.MR.v001_HNVCGBGXN_S3_318_HEK293T",
                     "VAI.MR.v001_HNVCGBGXN_S3_334_HEK293T", "VAI.MR.v001_HNVCGBGXN_S3_379_HEK293T")
p1_neg_controls %in% colnames(vasa)
# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

p2_neg_controls <- c("VAI.MR.v002_HNVCGBGXN_S4_009_K562", "VAI.MR.v002_HNVCGBGXN_S4_047_K562",
                     "VAI.MR.v002_HNVCGBGXN_S4_108_K562", "VAI.MR.v002_HNVCGBGXN_S4_123_K562",
                     "VAI.MR.v002_HNVCGBGXN_S4_207_HEK293T", "VAI.MR.v002_HNVCGBGXN_S4_318_HEK293T",
                     "VAI.MR.v002_HNVCGBGXN_S4_334_HEK293T", "VAI.MR.v002_HNVCGBGXN_S4_379_HEK293T")

p2_neg_controls %in% colnames(vasa)
# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

## grab just the HEK cells
vasa.hek <- vasa[,grep("HEK293T", colnames(vasa))]
vasa.hek <- cbind(vasa.hek, vasa[,c(677:680)])

# write.table(vasa.hek, file = "~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/vasa/te_quants/hek293t_vasa_clean_count_matrix_filtered_cpm_TE.tsv",
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## grab just the K562 cells
vasa.k562 <- vasa[,grep("K562", colnames(vasa))]
vasa.k562 <- cbind(vasa.k562, vasa[,c(677:680)])
# write.table(vasa.k562, file = "~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/vasa/te_quants/k562_vasa_clean_count_matrix_filtered_cpm_TE.tsv",
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

#dataSTORM.filt <- dataSTORM.filt[rowMeans(dataSTORM.filt[,c(1:104)]) > 1,]
# K562
#STORM_BulkTotal <- merge(BulkTotal, dataSTORM.filt[,c(1:105)], by = "Geneid")
#STORM_BulkTotal <- merge(BulkTotal, dataSTORM.filt[,c(1:104, 107)], by = "Geneid")

# HEK293T
# VASA
#STORM_BulkTotal <- merge(BulkTotal, dataSTORM.filt[,c(1:132)], by = "Geneid")
STORM_BulkTotal <- merge(BulkTotal, dataSTORM.filt[,c(1:132)], by = "Geneid")

## add in ss2
# STORM_BulkTotal_ss2 <- merge(STORM_BulkTotal, ss2[,c(1:64)], by = "Geneid")

## add in vasa
#STORM_BulkTotal_vasa <- merge(STORM_BulkTotal, vasa.k562[,c(1:326)], by = "Geneid")
STORM_BulkTotal_vasa <- merge(STORM_BulkTotal, vasa.hek[,c(1:352)], by = "Geneid")

## add in sstotal
# STORM_BulkTotal_sstotal <- merge(STORM_BulkTotal, sstotal[,c(1:249)], by = "Geneid")

write.table(STORM_BulkTotal_vasa,
            file = "~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/bulk_storm_vasa_hek_merged_te_quants_cpm.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

## filter a bit
# cat bulk_storm_vasa_hek_merged_te_quants_cpm.txt | tail -n+2 | cut -f2,3,7-488 | awk '{c=0;for(i=1;i<=NF;++i){c+=$i};print c}' > rowsums.txt
full.filt <- read.delim("~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/rowsums.txt", header = FALSE)
STORM_BulkTotal_vasa.filt <- STORM_BulkTotal_vasa[full.filt$V1 > 0,]
## mega hek
STORM_BulkTotal_vasa_sstotal <- merge(STORM_BulkTotal_vasa.filt, sstotal[,c(1:249)], by = "Geneid")

write.table(STORM_BulkTotal_vasa_sstotal,
            file = "~/Documents/manuscripts/mini_scrna_2020/te_vasa_sstotal/bulk_storm_vasa_sstotal_hek_merged_te_quants_cpm.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

# all single cells
df.OG2 <- data.matrix(STORM_BulkTotal_vasa_sstotal[STORM_BulkTotal_vasa_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:736)])
colnames(df.OG2) <- colnames(STORM_BulkTotal_vasa_sstotal[STORM_BulkTotal_vasa_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:736)])
rownames(df.OG2) <- STORM_BulkTotal_vasa_sstotal[STORM_BulkTotal_vasa_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,6,7:736)]$repFamily

# K562 SS2
# df.OG2 <- data.matrix(STORM_BulkTotal_ss2[STORM_BulkTotal_ss2$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:173)])
# colnames(df.OG2) <- colnames(STORM_BulkTotal_ss2[STORM_BulkTotal_ss2$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:173)])
# rownames(df.OG2) <- STORM_BulkTotal_ss2[STORM_BulkTotal_ss2$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,6,7:173)]$repFamily

#K562 VASA
# df.OG2 <- data.matrix(STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:435)])
# colnames(df.OG2) <- colnames(STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:435)])
# rownames(df.OG2) <- STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,6,7:435)]$repFamily

#HEK VASA
# df.OG2 <- data.matrix(STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:488)])
# colnames(df.OG2) <- colnames(STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:488)])
# rownames(df.OG2) <- STORM_BulkTotal_vasa[STORM_BulkTotal_vasa$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,6,7:488)]$repFamily

#HEK SSTOTAL
# df.OG2 <- data.matrix(STORM_BulkTotal_sstotal[STORM_BulkTotal_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:385)])
# colnames(df.OG2) <- colnames(STORM_BulkTotal_sstotal[STORM_BulkTotal_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,7:385)])
# rownames(df.OG2) <- STORM_BulkTotal_sstotal[STORM_BulkTotal_sstotal$repClass %in% c("LTR", "LINE", "SINE"), c(2,3,6,7:385)]$repFamily

library(ComplexHeatmap)
Heatmap(df.OG2, column_title = "Cells",name= "LTR, LINE, and SINE CPM",col = colorRamp2(c(0,1,2,3,4,5,6,7,8,9,10), c(rev(mypalette5))), 
        cluster_rows = T, cluster_columns = T,show_row_names = F, row_split = factor(rownames(df.OG2)))

## not bad

# bulk_mean_te_exp <- rowMeans(log1p(df.OG2.erv.line[,c(1:2)]))
# sc_mean_te_exp <- rowMeans(log1p(df.OG2.erv.line[,c(3:106)]))

## grab loci that are covered by at least 50% of cells
# df.OG2.erv.line.filt <- df.OG2.erv.line[rowSums(apply(df.OG2.erv.line[,c(1:106)] > 0, 2, as.numeric))/106 > 0.5,]
## K562
#df.OG2.filt <- df.OG2[rowSums(apply(df.OG2[,c(1:106)] > 0, 2, as.numeric))/106 > 0.5,]
## HEK293T
## find loci covered in 50% of cells, sampling the VASA data to similar numbers as STORM

# K562
## sample 63 STORM cells
set.seed(1988)
storm.samp <- sample(colnames(df.OG2)[3:106], size = 63, replace = FALSE)
ss2.cells <- colnames(df.OG2)[107:169]
bulk.samps <- colnames(df.OG2)[1:2]

# K562
## sample 104 VASA cells
# set.seed(1988)
# vasa.samp <- sample(colnames(df.OG2)[107:431], size = 104, replace = FALSE)
# storm.cells <- colnames(df.OG2)[3:106]
# bulk.samps <- colnames(df.OG2)[1:2]

# HEK
## sample 130 VASA cells
# set.seed(1988)
# vasa.samp <- sample(colnames(df.OG2)[134:484], size = 130, replace = FALSE)
# storm.cells <- colnames(df.OG2)[3:133]
# bulk.samps <- colnames(df.OG2)[1:2]

# HEK
## sample 131 SSTOTAL cells
# set.seed(1988)
# sstotal.samp <- sample(colnames(df.OG2)[134:381], size = 131, replace = FALSE)
# storm.cells <- colnames(df.OG2)[3:133]
# bulk.samps <- colnames(df.OG2)[1:2]

## label tes in the rownames
rownames(df.OG2) <- paste0("te_id_", 1:nrow(df.OG2), "_", rownames(df.OG2))

df.OG2.sub <- df.OG2[,colnames(df.OG2) %in% c(bulk.samps,
                                              storm.samp,
                                              ss2.cells)]

# df.OG2.sub <- df.OG2[,colnames(df.OG2) %in% c(bulk.samps,
#                                               storm.cells,
#                                               vasa.samp)]
# df.OG2.sub <- df.OG2[,colnames(df.OG2) %in% c(bulk.samps,
#                                               storm.cells,
#                                               sstotal.samp)]
df.OG2.sub.filt <- df.OG2.sub[rowSums(apply(df.OG2.sub[,c(1:128)] > 0, 2, as.numeric))/128 >= 0.25,]

df.OG2.filt <- df.OG2[rownames(df.OG2) %in% rownames(df.OG2.sub.filt),]


## all cells
df.OG2.filt <- df.OG2[rowSums(apply(df.OG2[,c(1:732)] > 0, 2, as.numeric))/732 >= 0.1,]

## filter out some low quality cells
## K562 STORM remove low qual cells
# df.OG2.filt <- df.OG2.filt[,!colnames(df.OG2.filt) %in% c("K4", "D12", "P6", "D11")]


## all HEK cells
bulk_samps <- c("ERR304487", "ERR304488")
storm_cells <- colnames(df.OG2.filt)[3:133]
vasa_cells <- colnames(df.OG2.filt)[134:484]
sstotal_cells <- colnames(df.OG2.filt)[485:732]

## K562 SS2
# ss2_cell_ids <- colnames(ss2)[1:63]

## K562 Vasa
# vasa_cell_ids <- colnames(df.OG2.filt)[107:431]

## HEK293T Vasa
# vasa_cell_ids <- colnames(df.OG2.filt)[134:484]

## HEK293T SSTotal
# sstotal_cell_ids <- colnames(df.OG2.filt)[134:381]

## K562
# col_split <- as.factor(ifelse(colnames(df.OG2.filt) %in% c("K5621", "K5622"),
#                               "Bulk",
#                               ifelse(colnames(df.OG2.filt) %in% ss2_cell_ids,
#                                      "SMART-seq2",
#                                      "STORM-seq")))
# VASA K562
# col_split <- as.factor(ifelse(colnames(df.OG2.filt) %in% c("K5621", "K5622"),
#                               "Bulk",
#                               ifelse(colnames(df.OG2.filt) %in% vasa_cell_ids,
#                                      "VASA-seq",
#                                      "STORM-seq")))


# all cells HEK293T
col_split <- as.factor(ifelse(colnames(df.OG2.filt) %in% bulk_samps,
                              "Bulk",
                              ifelse(colnames(df.OG2.filt) %in% vasa_cells,
                                     "VASA-seq",
                                     ifelse(colnames(df.OG2.filt) %in% storm_cells,
                                            "STORM-seq",
                                            "SMART-seq Total"))))

# VASA HEK293T
# col_split <- as.factor(ifelse(colnames(df.OG2.filt) %in% c("ERR304487", "ERR304488"),
#                               "Bulk",
#                               ifelse(colnames(df.OG2.filt) %in% vasa_cell_ids,
#                                      "VASA-seq",
#                                      "STORM-seq")))

# SSTotal HEK293T
# col_split <- as.factor(ifelse(colnames(df.OG2.filt) %in% c("ERR304487", "ERR304488"),
#                               "Bulk",
#                               ifelse(colnames(df.OG2.filt) %in% sstotal_cell_ids,
#                                      "SMART-seq Total",
#                                      "STORM-seq")))
#   
## K562 with SS2
# Heatmap(log1p(df.OG2.filt), column_title = "K-562",name= "LTR, LINE, and SINE CPM",col = colorRamp2(c(0,1,2,3,4,5,6,7,8,9,10), c(rev(mypalette5))), 
#         cluster_rows = T, cluster_columns = T,show_row_names = F, row_split = factor(rownames(df.OG2.filt)),
#         column_split = col_split, show_column_names = F)

## HEK293T with Vasa
# Heatmap(log1p(df.OG2.filt), column_title = "HEK293T",name= "LTR, LINE, and SINE CPM",col = colorRamp2(c(0,1,2,3,4,5,6,7,8,9,10), c(rev(mypalette5))), 
#         cluster_rows = T, cluster_columns = T,show_row_names = F, row_split = factor(rownames(df.OG2.filt)),
#         column_split = col_split, show_column_names = F)

  
# looks decent

## all cells
bulk_mean_te_exp <- rowMeans(df.OG2.filt[,1:2])
storm_mean_te_exp <- rowMeans(df.OG2.filt[,c(3:133)])
vasa_mean_te_exp <- rowMeans(df.OG2.filt[,c(134:484)])
sstotal_mean_te_exp <- rowMeans(df.OG2.filt[,c(485:732)])

## what's the correlations with bulk?
# bulk_mean_te_exp <- rowMeans(df.OG2.filt[,1:2])
# storm_mean_te_exp <- rowMeans(df.OG2.filt[,c(3:106)])
# ss2_mean_te_exp <- rowMeans(df.OG2.filt[,c(107:169)])

# K562 VASA
# bulk_mean_te_exp <- rowMeans(df.OG2.filt[,1:2])
# storm_mean_te_exp <- rowMeans(df.OG2.filt[,c(3:106)])
# vasa_mean_te_exp <- rowMeans(df.OG2.filt[,c(107:431)])


# HEK VASA
# bulk_mean_te_exp <- rowMeans(df.OG2.filt[,1:2])
# storm_mean_te_exp <- rowMeans(df.OG2.filt[,c(3:133)])
# vasa_mean_te_exp <- rowMeans(df.OG2.filt[,c(134:484)])

# HEK SSTotal
# bulk_mean_te_exp <- rowMeans(df.OG2.filt[,1:2])
# storm_mean_te_exp <- rowMeans(df.OG2.filt[,c(3:133)])
# sstotal_mean_te_exp <- rowMeans(df.OG2.filt[,c(134:381)])

plot(log1p(bulk_mean_te_exp), log1p(storm_mean_te_exp))
# plot(log1p(bulk_mean_te_exp), log1p(ss2_mean_te_exp))
plot(log1p(bulk_mean_te_exp), log1p(vasa_mean_te_exp))
plot(log1p(bulk_mean_te_exp), log1p(sstotal_mean_te_exp))


to_plot <- data.frame(Bulk = log1p(bulk_mean_te_exp),
                      STORM = log1p(storm_mean_te_exp),
                      `VASA-seq` = log1p(vasa_mean_te_exp),
                      `SMART-seq Total` = log1p(sstotal_mean_te_exp))

# to_plot <- data.frame(Bulk = log1p(bulk_mean_te_exp),
#                       STORM = log1p(storm_mean_te_exp),
#                       `SMART-seq2` = log1p(ss2_mean_te_exp))

# to_plot <- data.frame(Bulk = log1p(bulk_mean_te_exp),
#                       STORM = log1p(storm_mean_te_exp),
#                       `VASA-seq` = log1p(vasa_mean_te_exp))

# to_plot <- data.frame(Bulk = log1p(bulk_mean_te_exp),
#                       STORM = log1p(storm_mean_te_exp),
#                       `SMART-seq Total` = log1p(sstotal_mean_te_exp))

to_plot.storm_vs_bulk <- to_plot[,c("Bulk", "STORM")]
to_plot.storm_vs_bulk$group <- "Bulk_vs_STORM"
# to_plot.ss2_vs_bulk <- to_plot[,c("Bulk", "SMART.seq2")]
# to_plot.ss2_vs_bulk$group <- "Bulk_vs_SMART.seq2"

to_plot.vasa_vs_bulk <- to_plot[,c("Bulk", "VASA.seq")]
to_plot.vasa_vs_bulk$group <- "Bulk_vs_VASA.seq"

to_plot.sstotal_vs_bulk <- to_plot[,c("Bulk", "SMART.seq.Total")]
to_plot.sstotal_vs_bulk$group <- "Bulk_vs_SMART.seq.Total"

colnames(to_plot.storm_vs_bulk)[2] <- "SingleCell"
#colnames(to_plot.ss2_vs_bulk)[2] <- "SingleCell"
colnames(to_plot.vasa_vs_bulk)[2] <- "SingleCell"
colnames(to_plot.sstotal_vs_bulk)[2] <- "SingleCell"


to_plot_all <- rbind(to_plot.storm_vs_bulk,
                     to_plot.vasa_vs_bulk,
                     to_plot.sstotal_vs_bulk)

# to_plot_all <- rbind(to_plot.storm_vs_bulk,
#                      to_plot.ss2_vs_bulk)

# to_plot_all <- rbind(to_plot.storm_vs_bulk,
#                      to_plot.vasa_vs_bulk)

# to_plot_all <- rbind(to_plot.storm_vs_bulk,
#                      to_plot.sstotal_vs_bulk)

#make some nicer plots for correlation

# All cells
cor(to_plot$Bulk, to_plot$STORM)
# 0.8094185
cor(to_plot$Bulk, to_plot$VASA.seq)
# 0.4029486
cor(to_plot$Bulk, to_plot$SMART.seq.Total)
# 0.06020767

# cor(to_plot$Bulk, to_plot$STORM)
# # 0.7622346
# cor(to_plot$Bulk, to_plot$SMART.seq2)
# # 0.2678593

# VASA K562
# cor(to_plot$Bulk, to_plot$STORM)
# # 0.8088748
# cor(to_plot$Bulk, to_plot$VASA.seq)
# # 0.1631234

# VASA HEK
# cor(to_plot$Bulk, to_plot$STORM)
# 0.8063334
# cor(to_plot$Bulk, to_plot$VASA.seq)
# -0.01955761

# SSTotal HEK
# cor(to_plot$Bulk, to_plot$STORM)
# 0.8750064
# cor(to_plot$Bulk, to_plot$SMART.seq.Total)
# 0.02974612

ggplot(to_plot_all, aes(x = Bulk, y = SingleCell, color = group)) +
  geom_point(alpha = 0.45, size = 2) +
  geom_smooth(method = lm) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(`Bulk_vs_STORM` = "#63197F",
                                  `Bulk_vs_SMART.seq2` = "#DE4766"),
                          labels = c("SMART-seq2",
                                     "STORM-seq")) +
  scale_x_continuous(breaks = c(0,3,6,9),
                     limits = c(0,9)) +
  scale_y_continuous(breaks = c(0,3,6,9),
                     limits = c(0,9)) +
  ylab("Single-cell TE expression\nlog(CPM+1)") +
  xlab("Bulk TE expression\nlog(CPM+1)") +
  theme_bw(12) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

# ggplot(to_plot_all, aes(x = Bulk, y = SingleCell, color = group)) +
#   geom_point(alpha = 0.45, size = 2) +
#   geom_smooth(method = lm) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c(`Bulk_vs_STORM` = "#63197F",
#                                 `Bulk_vs_VASA.seq` = "#CF5917"),
#                         labels = c("STORM-seq",
#                                    "VASA-seq")) +
#   scale_x_continuous(breaks = c(0,3,6,9),
#                      limits = c(0,9)) +
#   scale_y_continuous(breaks = c(0,3,6,9),
#                      limits = c(0,9)) +
#   ylab("Single-cell TE expression\nlog(CPM+1)") +
#   xlab("Bulk TE expression\nlog(CPM+1)") +
#   theme_bw(12) +
#   theme(
#     legend.title = element_blank(),
#     legend.text = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16)
#   )

# ggplot(to_plot_all, aes(x = Bulk, y = SingleCell, color = group)) +
#   geom_point(alpha = 0.45, size = 2) +
#   geom_smooth(method = lm) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c(`Bulk_vs_STORM` = "#63197F",
#                                 `Bulk_vs_SMART.seq.Total` = "#8FC2FD"),
#                      labels = c("SMART-seq Total",
#                                 "STORM-seq")) +
#   scale_x_continuous(breaks = c(0,3,6,9),
#                      limits = c(0,9)) +
#   scale_y_continuous(breaks = c(0,3,6,9),
#                      limits = c(0,9)) +
#   ylab("Single-cell TE expression\nlog(CPM+1)") +
#   xlab("Bulk TE expression\nlog(CPM+1)") +
#   theme_bw(12) +
#   theme(
#     legend.title = element_blank(),
#     legend.text = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16)
#   )

ggplot(to_plot_all, aes(x = Bulk, y = SingleCell, color = group)) +
  geom_point(alpha = 0.45, size = 2) +
  geom_smooth(method = lm) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(`Bulk_vs_STORM` = "#63197F",
                                `Bulk_vs_VASA.seq` = "#CF5917",
                                `Bulk_vs_SMART.seq.Total` = "#8FC2FD"),
                        labels = c("SMART-seq Total",
                                   "STORM-seq",
                                   "VASA-seq")) +
  scale_x_continuous(breaks = c(0,3,6,9),
                     limits = c(0,9)) +
  scale_y_continuous(breaks = c(0,3,6,9),
                     limits = c(0,9)) +
  ylab("Single-cell TE expression\nlog(CPM+1)") +
  xlab("Bulk TE expression\nlog(CPM+1)") +
  facet_wrap(~group) +
  theme_bw(12) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
