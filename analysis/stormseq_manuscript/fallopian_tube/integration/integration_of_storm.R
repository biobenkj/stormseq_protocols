## Integrate the two patients' STORM-seq data
library(scran)
library(scater)
library(scuttle)
library(BiocSingular)
library(cowplot)

## import data
storm_pat_1 <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")
storm_pat_2 <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")

## do we need to do the usual re-scaling with Batchelor?
library(batchelor)
## rescale
rescaled_data <- multiBatchNorm(storm_pat_1,
                                storm_pat_2)
rescale_storm_pat_1 <- rescaled_data[[1]]
rescale_storm_pat_2 <- rescaled_data[[2]]

## combine variance
storm_pat_1.dec <- modelGeneVar(storm_pat_1)
storm_pat_2.dec <- modelGeneVar(storm_pat_2)

combined_var <- combineVar(storm_pat_1.dec,
                           storm_pat_2.dec)

chosen.hvgs <- combined_var$bio > 0
sum(chosen.hvgs)
# 17573

top.hvgs <- scran::getTopHVGs(combined_var, prop=0.2)

## plot out the uncorrected batches
## need to toss some metadata in colData to cbind
colData(rescale_storm_pat_1) <- colData(rescale_storm_pat_1)[,c(1:4)]
colData(rescale_storm_pat_2) <- colData(rescale_storm_pat_2)[,c(1:4)]

rescale_storm_pat_1$batch <- "NovaSeq"
rescale_storm_pat_2$batch <- "HiSeq"

## throw out some reducedDims too
reducedDim(rescale_storm_pat_1, "PCAsub") <- NULL
reducedDim(rescale_storm_pat_2, "PCAsub") <- NULL

reducedDim(rescale_storm_pat_1, "denSNE") <- NULL
reducedDim(rescale_storm_pat_2, "denSNE") <- NULL

reducedDim(rescale_storm_pat_1, "densMAP") <- NULL
reducedDim(rescale_storm_pat_2, "densMAP") <- NULL

rescale_storm_all_pats <- cbind(rescale_storm_pat_1,
                                rescale_storm_pat_2)
## runPCA on uncorrected values
rescale_storm_all_pats <- runPCA(rescale_storm_all_pats, subset_row=top.hvgs,
                                 BSPARAM=BiocSingular::IrlbaParam())

snn.gr <- buildSNNGraph(rescale_storm_all_pats, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=rescale_storm_all_pats$batch)
tab
rescale_storm_all_pats <- runTSNE(rescale_storm_all_pats, dimred="PCA")
rescale_storm_all_pats <- runUMAP(rescale_storm_all_pats, dimred="PCA")
plotTSNE(rescale_storm_all_pats, colour_by="batch")
plotUMAP(rescale_storm_all_pats, colour_by="batch")

## eh, maybe... still do MNN and see
set.seed(1988)
mnn.out <- fastMNN(rescale_storm_pat_1, rescale_storm_pat_2, d = 50, k=15, subset.row=chosen.hvgs,
                   BSPARAM=BiocSingular::IrlbaParam())
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn
mnn.out$cluster <- clusters.mnn
## yeah that's better!
library(densvis)
mnn.out.densmap <- densmap(reducedDim(mnn.out, "corrected"), n_components = 3L,
                           n_neighbors = 15L, metric = "euclidean")
#mnn.out <- runTSNE(mnn.out, dimred = "corrected")
reducedDim(mnn.out, "densMAP") <- mnn.out.densmap
plotReducedDim(mnn.out, colour_by = "cluster", dimred = "densMAP", ncomponents = 2)
plotTSNE(mnn.out, colour_by = "batch")

to_plot <- data.frame(Dim_1 = reducedDim(mnn.out, "TSNE")[,1],
                      Dim_2 = reducedDim(mnn.out, "TSNE")[,2],
                      Patient = c(rep("Patient 1", 308),
                                  rep("Patient 2", 333)))
## tack back on the TPM to show cell types are retained in the joint embedding
to_plot$epcam <- c(storm_pat_1$epcam, storm_pat_2$epcam)
colData(storm_pat_1)$tmem173 <- assays(storm_pat_1)$izar[rowData(storm_pat_1)$symbol == "STING1",]
colData(storm_pat_2)$tmem173 <- assays(storm_pat_2)$izar[rowData(storm_pat_2)$symbol == "STING1",]
to_plot$sting1 <- c(storm_pat_1$tmem173, storm_pat_2$tmem173)
to_plot$secretory <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(storm_pat_2)[,c("pax8", "krt7")])))
to_plot$ciliated <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("foxj1", "ccdc17", "ccdc78")])),
                       rowSums(as.data.frame(colData(storm_pat_2)[,c("foxj1", "ccdc17", "ccdc78")])))
to_plot$immune <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("cd45", "cd11c", "cd14")])),
                      rowSums(as.data.frame(colData(storm_pat_2)[,c("cd45", "cd11c", "cd14")])))
to_plot$cluster <- mnn.out$cluster
integrated <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = Patient)) +
  geom_point() +
  scale_color_manual(values = c("maroon", "gray"), name = "Patient") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Integrated Patient FTE STORM-seq") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

epcam <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = epcam)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EpCAM gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

sting1 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = sting1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("STING1 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

secretory <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = secretory)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Secretory cell marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ciliated <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = ciliated)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Ciliated cell marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

immune <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = immune)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Immune marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )
plot_grid(integrated, epcam, secretory, ciliated, immune, ncol = 2, nrow = 3)

storm_mnn_clusts <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = factor(cluster))) +
  geom_point() +
  scale_color_viridis_d(name = "Cluster") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Integrated Patient FTE STORM-seq") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
storm_mnn_clusts

## integrate with SS2 data
yau_all_data <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/ss2_salmon_requant/ss2_salmon_ens102_requant_grch38_velocessor_normd_data_with_geneids_with_embeddings.rds")
names(yau_all_data) <- paste0("ss2_patient_", 1:4)

## pull some key gene expression for the Yau data
yau_all_data <- lapply(yau_all_data, function(x) {
  x$epcam <- assays(x)$izar[rowData(x)$symbol == "EPCAM",]
  x$foxj1 <- assays(x)$izar[rowData(x)$symbol == "FOXJ1",]
  x$ccdc17 <- assays(x)$izar[rowData(x)$symbol == "CCDC17",]
  x$ccdc78 <- assays(x)$izar[rowData(x)$symbol == "CCDC78",]
  x$cd45 <- assays(x)$izar[rowData(x)$symbol == "PTPRC",]
  x$cd11c <- assays(x)$izar[rowData(x)$symbol == "ITGAX",]
  x$cd14 <- assays(x)$izar[rowData(x)$symbol == "CD14",]
  x$krt7 <- assays(x)$izar[rowData(x)$symbol == "KRT7",]
  x$pax8 <- assays(x)$izar[rowData(x)$symbol == "PAX8",]
  x$cd34 <- assays(x)$izar[rowData(x)$symbol == "CD34",]
  return(x)
})

## combine the Yau data into 1 SCE
all_yau_data <- do.call(cbind, lapply(yau_all_data, function(x) {
  reducedDim(x, "PCAsub") <- NULL
  return(x)
}))



## rescale
rescaled_data <- multiBatchNorm(storm_pat_1,
                                storm_pat_2,
                                all_yau_data)

rescale_storm_pat_1 <- rescaled_data[[1]]
rescale_storm_pat_2 <- rescaled_data[[2]]

rescale_yau_all <- rescaled_data[[3]]

## combine variance
storm_pat_1.dec <- modelGeneVar(storm_pat_1)
storm_pat_2.dec <- modelGeneVar(storm_pat_2)

ss2_yau_all.dec <- modelGeneVar(all_yau_data)

combined_var <- combineVar(storm_pat_1.dec,
                           storm_pat_2.dec,
                           ss2_yau_all.dec)

chosen.hvgs <- combined_var$bio > 0
sum(chosen.hvgs)
# 16598

#top.hvgs <- scran::getTopHVGs(combined_var, prop=0.2)

## need to correct colData to cbind
colData(rescale_storm_pat_1) <- colData(rescale_storm_pat_1)[,c(1:4)]
colData(rescale_storm_pat_2) <- colData(rescale_storm_pat_2)[,c(1:4)]
reducedDim(rescale_storm_pat_1, "PCAsub") <- NULL
reducedDim(rescale_storm_pat_2, "PCAsub") <- NULL
reducedDim(rescale_storm_pat_1, "denSNE") <- NULL
reducedDim(rescale_storm_pat_2, "denSNE") <- NULL
reducedDim(rescale_storm_pat_1, "densMAP") <- NULL
reducedDim(rescale_storm_pat_2, "densMAP") <- NULL
rescale_storm_pat_1$batch <- "STORM_patient1"
rescale_storm_pat_2$batch <- "STORM_patient2"

colData(rescale_yau_all) <- colData(rescale_yau_all)[,c(1:4)]
# colData(rescale_ss2_pat_2) <- colData(rescale_ss2_pat_2)[,c(1:4)]
# colData(rescale_ss2_pat_3) <- colData(rescale_ss2_pat_3)[,c(1:4)]
# colData(rescale_ss2_pat_4) <- colData(rescale_ss2_pat_4)[,c(1:4)]
reducedDim(rescale_yau_all, "PCAsub") <- NULL
# reducedDim(rescale_ss2_pat_2, "PCAsub") <- NULL
# reducedDim(rescale_ss2_pat_3, "PCAsub") <- NULL
# reducedDim(rescale_ss2_pat_4, "PCAsub") <- NULL
reducedDim(rescale_yau_all, "denSNE") <- NULL
# reducedDim(rescale_ss2_pat_2, "denSNE") <- NULL
# reducedDim(rescale_ss2_pat_3, "denSNE") <- NULL
# reducedDim(rescale_ss2_pat_4, "denSNE") <- NULL
reducedDim(rescale_yau_all, "densMAP") <- NULL
# reducedDim(rescale_ss2_pat_2, "densMAP") <- NULL
# reducedDim(rescale_ss2_pat_3, "densMAP") <- NULL
# reducedDim(rescale_ss2_pat_4, "densMAP") <- NULL
rescale_yau_all$batch <- "SS2"
# rescale_ss2_pat_2$batch <- "SS2_patient2"
# rescale_ss2_pat_3$batch <- "SS2_patient3"
# rescale_ss2_pat_4$batch <- "SS2_patient4"

rescale_all_pats <- cbind(rescale_storm_pat_1,
                          rescale_storm_pat_2,
                          rescale_yau_all)

## runPCA on uncorrected values
rescale_all_pats <- runPCA(rescale_all_pats, subset_row=chosen.hvgs,
                                 BSPARAM=BiocSingular::IrlbaParam())

snn.gr <- buildSNNGraph(rescale_all_pats, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=rescale_all_pats$batch)
tab
rescale_all_pats <- runTSNE(rescale_all_pats, dimred="PCA")
rescale_all_pats <- runUMAP(rescale_all_pats, dimred="PCA")
plotTSNE(rescale_all_pats, colour_by="cluster")
plotUMAP(rescale_all_pats, colour_by="batch")

set.seed(1988)
mnn.out <- fastMNN(rescale_storm_pat_1, rescale_storm_pat_2,
                   rescale_yau_all,
                   d = 50, k = 15, subset.row=chosen.hvgs,
                   BSPARAM=BiocSingular::RandomParam(deferred = TRUE))
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn
mnn.out$cluster <- clusters.mnn
## yeah that's better!
library(densvis)
mnn.out.densmap <- densmap(reducedDim(mnn.out, "corrected"), n_components = 2L,
                           n_neighbors = 15L, metric = "euclidean")
mnn.out <- runTSNE(mnn.out, dimred = "corrected")
#mnn.out <- runUMAP(mnn.out, dimred = "corrected")
reducedDim(mnn.out, "densMAP") <- mnn.out.densmap
plotReducedDim(mnn.out, colour_by = "cluster", dimred = "densMAP", ncomponents = 2)
## go with the TSNE here
#mnn.out$batch <- rescale_all_pats$batch
plotTSNE(mnn.out, colour_by = "cluster")
#plotUMAP(mnn.out, colour_by = "batch")


mnn.out$batch.collapse <- c(rep("STORM-seq", 641),
                            rep("SMART-seq2", 1762))
to_plot <- data.frame(Dim_1 = reducedDim(mnn.out, "TSNE")[,1],
                      Dim_2 = reducedDim(mnn.out, "TSNE")[,2],
                      Patient = mnn.out$batch,
                      Platform = mnn.out$batch.collapse)
## tack back on the TPM to show cell types are retained in the joint embedding
to_plot$epcam <- c(storm_pat_1$epcam, storm_pat_2$epcam,
                   yau_all_data[[1]]$epcam, yau_all_data[[2]]$epcam,
                   yau_all_data[[3]]$epcam, yau_all_data[[4]]$epcam)

to_plot$secretory <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(storm_pat_2)[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(yau_all_data[[1]])[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(yau_all_data[[2]])[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(yau_all_data[[3]])[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(yau_all_data[[4]])[,c("pax8", "krt7")])))

to_plot$ciliated <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(storm_pat_2)[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(yau_all_data[[1]])[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(yau_all_data[[2]])[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(yau_all_data[[3]])[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(yau_all_data[[4]])[,c("foxj1", "ccdc17", "ccdc78")])))

to_plot$immune <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(storm_pat_2)[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(yau_all_data[[1]])[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(yau_all_data[[2]])[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(yau_all_data[[3]])[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(yau_all_data[[4]])[,c("cd45", "cd11c", "cd14")])))

library(gghighlight)
integrated <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = Platform)) +
  geom_point() +
  gghighlight(Platform == "STORM-seq") +
  scale_color_manual(values = c("maroon", "gray"), name = "Platform") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Patient 1 and 2 FTE STORM-seq") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

integrated2 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = Platform)) +
  geom_point() +
  gghighlight(Platform == "SMART-seq2") +
  scale_color_manual(values = c("black", "gray"), name = "Platform") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Patients 1-4 FTE SMART-seq2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

epcam <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = epcam)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EpCAM gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

immune <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = immune)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Immune marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

secretory <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = secretory)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Secretory cell marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ciliated <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = ciliated)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Ciliated cell marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(integrated, integrated2, epcam, secretory, ciliated, immune,
          ncol = 2, nrow = 3)


### marker gene discovery of integrated STORM patients 1 and 2
## collapse clusters
library(dplyr)

reducedDim(storm_pat_1, "PCAsub") <- NULL
reducedDim(storm_pat_2, "PCAsub") <- NULL

reducedDim(storm_pat_1, "denSNE") <- NULL
reducedDim(storm_pat_2, "denSNE") <- NULL

reducedDim(storm_pat_1, "densMAP") <- NULL
reducedDim(storm_pat_2, "densMAP") <- NULL

storm_pat_1$batch <- "patient1"
storm_pat_2$batch <- "patient2"
storm_all_pats <- cbind(storm_pat_1,
                        storm_pat_2)

markers.sce.cluster.nocollapse <- scran::findMarkers(storm_all_pats, mnn.out$cluster, block = storm_all_pats$batch,
                                                     test = "wilcox", direction = "up",
                                                     gene.names = rowData(storm_all_pats)$symbol,
                                                     pval.type = "some", assay.type = "izar")

## write out tables
lapply(names(markers.sce.cluster.nocollapse), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  write.table(markers,
              file = paste0("~/Documents/manuscripts/metabolic_flow_fte_2020/marker_expression_tables/cluster_",
                            x,
                            "_marker_exp_table.txt"),
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE,
              sep = '\t')
})

markers.sce.cluster.nocollapse2 <- scran::findMarkers(storm_pat_1, storm_pat_1$cluster,
                                                     test = "binom",
                                                     gene.names = rowData(storm_pat_1)$symbol,
                                                     pval.type = "some", assay.type = "izar")


marker_genes_by_cluster.nocollapse <- do.call(rbind, lapply(unique(mnn.out$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.05,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:11]))
  #auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  #markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

top_marker_genes_by_cluster.nocollapse <- do.call(rbind, lapply(unique(mnn.out$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.05,]
  ## subset to upper quartile of AUC
  #markers$AUC.ave <- rowMeans(as.matrix(markers[,4:11]))
  #auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  #markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- markers[order(markers$summary.AUC, decreasing = T),]
  if (nrow(markers) >= 20) {
    markers <- markers[1:20,]
  }
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

top_marker_genes_by_cluster.nocollapse2 <- do.call(rbind, lapply(unique(mnn.out$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse2[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.05,]
  markers <- markers[order(markers$summary.logFC, decreasing = T),]
  if (nrow(markers) >= 10) {
    markers <- markers[1:10,]
  }
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

top_marker_genes_by_cluster.nocollapse.single <- do.call(rbind, lapply(unique(mnn.out$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.05,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:11]))
  auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- markers[order(markers$AUC.ave, decreasing = T),]
  markers <- markers[1:5,]
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

## rename clusters
clust_names_nocollapse <- mnn.out$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "Secretory_arm",
    . == 2 ~ "Ciliated_arm",
    . == 3 ~ "Root_like",
    . == 4 ~ "Ciliated_1",
    . == 5 ~ "Secretory_1",
    . == 6 ~ "Ciliated_2",
    . == 7 ~ "Secretory_2",
    . == 8 ~ "Intermediate",
    . == 9 ~ "Immune"
  ))

storm_all_pats$cluster.names <- clust_names_nocollapse$cluster_cat

# dedupe
top_marker_genes_by_cluster.nocollapse <- subset(top_marker_genes_by_cluster.nocollapse, !duplicated(top_marker_genes_by_cluster.nocollapse$marker_genes))

## make a heatmap
storm_all_pats.markers.top <- storm_all_pats[rowData(storm_all_pats)$symbol %in% top_marker_genes_by_cluster.nocollapse$marker_genes,]
library(ComplexHeatmap)

rownames(storm_all_pats.markers.top) <- rowData(storm_all_pats.markers.top)$symbol

storm_all_pats.markers.top <- storm_all_pats.markers.top[grep("SNOR", rownames(storm_all_pats.markers.top),
                                                              invert = TRUE),]

Heatmap(assay(storm_all_pats.markers.top, "izar"),
        name = "log2(TPM+1)",
        column_split = storm_all_pats$cluster.names,
        top_annotation = HeatmapAnnotation(Patient = storm_all_pats.markers.top$batch, col = list(Patient = c(patient1 = "lightgray", patient2 = "darkgray"))),
        show_column_names = FALSE,
        show_row_names = TRUE, col = getJetColors(circlize.cols = TRUE),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))

## make some violin plots
## known EMT
# GATA2

storm_all_pats2 <- storm_all_pats
rownames(storm_all_pats2) <- rowData(storm_all_pats2)$symbol

colData(storm_all_pats2)$timp3 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "TIMP3",]
colData(storm_all_pats2)$dcn <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "DCN",]
colData(storm_all_pats2)$sox7 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "SOX7",]
colData(storm_all_pats2)$gata2 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "GATA2",]
colData(storm_all_pats2)$cavin2 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "CAVIN2",]
colData(storm_all_pats2)$zeb1 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "ZEB1",]
colData(storm_all_pats2)$sparc <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "SPARC",]
colData(storm_all_pats2)$tet1 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "TET1",]
colData(storm_all_pats2)$a2m <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "A2M",]
colData(storm_all_pats2)$il33 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "IL33",]
colData(storm_all_pats2)$tmem273 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "TMEM273",]
colData(storm_all_pats2)$bmp2 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "BMP2",]

## intermediate cells
colData(storm_all_pats2)$rims4 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "RIMS4",]
colData(storm_all_pats2)$adcy8 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "ADCY8",]
colData(storm_all_pats2)$catip <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "CATIP",]
colData(storm_all_pats2)$apobec4 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "APOBEC4",]
colData(storm_all_pats2)$snora73B <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "SNORA73B",]

## secretory cells
colData(storm_all_pats2)$crisp3 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "CRISP3",]
colData(storm_all_pats2)$msln <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "MSLN",]
colData(storm_all_pats2)$rspo1 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "RSPO1",]

## ciliated cells
colData(storm_all_pats2)$dnah9 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "DNAH9",]
colData(storm_all_pats2)$cfap54 <- assays(storm_all_pats2)$izar[rowData(storm_all_pats2)$symbol == "CFAP54",]

vio_plot <- data.frame(cluster_id = as.factor(mnn.out$cluster),
                       cluster_rename = storm_all_pats2$cluster.names,
                       batch = mnn.out$batch,
                       epcam = storm_all_pats2$epcam,
                       runx3 = storm_all_pats2$runx3,
                       pecam1 = storm_all_pats2$pecam1,
                       ovgp1 = storm_all_pats2$ovgp1,
                       cd45 = storm_all_pats2$cd45,
                       acta2 = storm_all_pats2$acta2,
                       timp3= storm_all_pats2$timp3,
                       dcn = storm_all_pats2$dcn,
                       sox7 = storm_all_pats2$sox7,
                       gata2 = storm_all_pats2$gata2,
                       cavin2 = storm_all_pats2$cavin2,
                       krt7 = storm_all_pats2$krt7,
                       pax8 = storm_all_pats2$pax8,
                       zeb1 = storm_all_pats2$zeb1,
                       sparc = storm_all_pats2$sparc,
                       tet1 = storm_all_pats2$tet1,
                       a2m = storm_all_pats2$a2m,
                       il33 = storm_all_pats2$il33,
                       tmem273 = storm_all_pats2$tmem273,
                       bmp2 = storm_all_pats2$bmp2,
                       tubb4 = storm_all_pats2$tubb4,
                       cd44 = storm_all_pats2$cd44,
                       itga6 = storm_all_pats2$itga6,
                       peg = rowSums(cbind(storm_all_pats2$epcam,
                                           storm_all_pats2$cd44,
                                           storm_all_pats2$itga6)),
                       rims4 = storm_all_pats2$rims4,
                       adcy8 = storm_all_pats2$adcy8,
                       catip = storm_all_pats2$catip,
                       apobec4 = storm_all_pats2$apobec4,
                       snora73B = storm_all_pats2$snora73B,
                       crisp3 = storm_all_pats2$crisp3,
                       rspo1 = storm_all_pats2$rspo1,
                       msln = storm_all_pats2$msln,
                       dnah9 = storm_all_pats2$dnah9,
                       cfap54 = storm_all_pats2$cfap54)

## intermediate cells
adcy8.plt <- ggplot(vio_plot, aes(x = cluster_id, y = adcy8,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("ADCY8") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
rims4.plt <- ggplot(vio_plot, aes(x = cluster_id, y = rims4,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("RIMS4") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
apobec4.plt <- ggplot(vio_plot, aes(x = cluster_id, y = apobec4,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("APOBEC4") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
catip.plt <- ggplot(vio_plot, aes(x = cluster_id, y = catip,
                                    group = cluster_id,
                                    color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CATIP") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
plot_grid(catip.plt,
          adcy8.plt,
          nrow = 1,
          ncol = 2)

##secretory
crisp3.plt <- ggplot(vio_plot, aes(x = cluster_id, y = crisp3,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CRISP3") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
ovgp1.plt <- ggplot(vio_plot, aes(x = cluster_id, y = ovgp1,
                                   group = cluster_id,
                                   color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("OVGP1") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
krt7.plt <- ggplot(vio_plot, aes(x = cluster_id, y = krt7,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("KRT7") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
msln.plt <- ggplot(vio_plot, aes(x = cluster_id, y = msln,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("MSLN") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
rspo1.plt <- ggplot(vio_plot, aes(x = cluster_id, y = rspo1,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("RSPO1") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
plot_grid(msln.plt,
          rspo1.plt,
          nrow = 1,
          ncol = 2)

## ciliated cells
cfap54.plt <- ggplot(vio_plot, aes(x = cluster_id, y = cfap54,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CFAP54") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
dnah9.plt <- ggplot(vio_plot, aes(x = cluster_id, y = dnah9,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("DNAH9") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
plot_grid(dnah9.plt,
          cfap54.plt,
          nrow = 1,
          ncol = 2)

## immune
cd45.plt <- ggplot(vio_plot, aes(x = cluster_id, y = cd45,
                                   group = cluster_id,
                                   color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CD45") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
runx3.plt <- ggplot(vio_plot, aes(x = cluster_id, y = runx3,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("RUNX3") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
plot_grid(cd45.plt,
          runx3.plt,
          nrow = 1,
          ncol = 2)

## cluster 3 progenitor like cells
gata2.plt <- ggplot(vio_plot, aes(x = cluster_id, y = gata2,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("GATA2") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

tmem273.plt <- ggplot(vio_plot, aes(x = cluster_id, y = tmem273,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("TMEM273") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

bmp2.plt <- ggplot(vio_plot, aes(x = cluster_id, y = bmp2,
                                    group = cluster_id,
                                    color = cluster_id)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("BMP2") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

timp3.plt <- ggplot(vio_plot, aes(x = cluster_id, y = timp3,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("TIMP3") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

dcn.plt <- ggplot(vio_plot, aes(x = cluster_id, y = dcn,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("DCN") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

sparc.plt <- ggplot(vio_plot, aes(x = cluster_id, y = sparc,
                                group = cluster_id,
                                color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("SPARC") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

zeb1.plt <- ggplot(vio_plot, aes(x = cluster_id, y = zeb1,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("ZEB1") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

acta2.plt <- ggplot(vio_plot, aes(x = cluster_id, y = acta2,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("ACTA2") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

cavin2.plt <- ggplot(vio_plot, aes(x = cluster_id, y = cavin2,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CAVIN2") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

il33.plt <- ggplot(vio_plot, aes(x = cluster_id, y = il33,
                                   group = cluster_id,
                                   color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("IL33") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

sox7.plt <- ggplot(vio_plot, aes(x = cluster_id, y = sox7,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("SOX7") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

plot_grid(gata2.plt, timp3.plt, cavin2.plt, zeb1.plt, bmp2.plt, nrow = 1, ncol = 5)

## PEG cells
epcam.plt <- ggplot(vio_plot, aes(x = cluster_id, y = epcam,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("EpCAM") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
cd44.plt <- ggplot(vio_plot, aes(x = cluster_id, y = cd44,
                                  group = cluster_id,
                                  color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("CD44") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
itga6.plt <- ggplot(vio_plot, aes(x = cluster_id, y = itga6,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("ITGA6") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )
peg.plt <- ggplot(vio_plot, aes(x = cluster_id, y = peg,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("PEG (EpCAM/CD44/ITGA6)") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

plot_grid(epcam.plt, cd44.plt,
          itga6.plt, peg.plt,
          nrow = 2, ncol = 2)

krt7.plt <- ggplot(vio_plot, aes(x = cluster_id, y = krt7,
                                 group = cluster_id,
                                 color = cluster_id)) +
  geom_violin(fill = "gray90", scale = "width", width = 0.8) +
  geom_jitter(shape = 16) +
  scale_color_viridis_d() +
  theme_half_open() +
  ggtitle("KRT7") +
  xlab("Cluster") +
  ylab("Expression log2(TPM+1)") +
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5)
  )

## make some marker gene expression plots
to_plot <- data.frame(Dim_1 = reducedDim(mnn.out, "TSNE")[,1],
                      Dim_2 = reducedDim(mnn.out, "TSNE")[,2],
                      Patient = c(rep("Patient 1", 308),
                                  rep("Patient 2", 333)))
## tack back on the TPM to show cell types are retained in the joint embedding
to_plot$epcam <- c(storm_pat_1$epcam, storm_pat_2$epcam)

colData(storm_pat_1)$runx3 <- assays(storm_pat_1)$izar[rowData(storm_pat_1)$symbol == "RUNX3",]
colData(storm_pat_2)$runx3 <- assays(storm_pat_2)$izar[rowData(storm_pat_2)$symbol == "RUNX3",]
to_plot$runx3 <- c(storm_pat_1$runx3, storm_pat_2$runx3)

colData(storm_pat_1)$acta2 <- assays(storm_pat_1)$izar[rowData(storm_pat_1)$symbol == "ACTA2",]
colData(storm_pat_2)$acta2 <- assays(storm_pat_2)$izar[rowData(storm_pat_2)$symbol == "ACTA2",]
colData(storm_pat_1)$timp3 <- assays(storm_pat_1)$izar[rowData(storm_pat_1)$symbol == "TIMP3",]
colData(storm_pat_2)$timp3 <- assays(storm_pat_2)$izar[rowData(storm_pat_2)$symbol == "TIMP3",]
to_plot$emt <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("acta2", "timp3")])),
                       rowSums(as.data.frame(colData(storm_pat_2)[,c("acta2", "timp3")])))

colData(storm_pat_1)$tmem173 <- assays(storm_pat_1)$izar[rowData(storm_pat_1)$symbol == "STING1",]
colData(storm_pat_2)$tmem173 <- assays(storm_pat_2)$izar[rowData(storm_pat_2)$symbol == "STING1",]
to_plot$sting1 <- c(storm_pat_1$tmem173, storm_pat_2$tmem173)
to_plot$secretory <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("pax8", "krt7")])),
                       rowSums(as.data.frame(colData(storm_pat_2)[,c("pax8", "krt7")])))
to_plot$ciliated <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("foxj1", "ccdc17", "ccdc78")])),
                      rowSums(as.data.frame(colData(storm_pat_2)[,c("foxj1", "ccdc17", "ccdc78")])))
to_plot$immune <- c(rowSums(as.data.frame(colData(storm_pat_1)[,c("cd45", "cd11c", "cd14")])),
                    rowSums(as.data.frame(colData(storm_pat_2)[,c("cd45", "cd11c", "cd14")])))
to_plot$cluster <- mnn.out$cluster
runx3 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = runx3)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("RUNX3 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

immune <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = immune)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Immune cell gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(immune, runx3)

emt <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = emt)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2) +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EMT marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

## get some proportion info
table(storm_all_pats2$cluster.names)
# Ciliated_1    Ciliated_2  Ciliated_arm        Immune 
# 33            65            57            16 
# Intermediate     Root_like   Secretory_1   Secretory_2 
# 38            58           105           194 
# Secretory_arm 
# 75 

## total 641
## ciliated including intermediate arm: 0.2418097
## secretory including intermediate arm: 0.5834633
## root like: 0.09048362
## intermediate: 0.05928237
## immune: 0.024961

## look at the cancer cell marker genes for EMT and Krt17 clusters
emt <- c("ACTA2", "TIMP3", "A2M", "TAGLN", "TPM2", "SPARC",
         "DCN", "MFAP4", "CRISPLD2", "LGALS1", "MYH11", "SPARCL1",
         "SLC2A3", "UBB", "SFRP4")
krt17 <- c("KRT17", "SNCG", "PIGR", "PLAT", "GPX3", "CHI3L1",
           "HLA-DPA1", "TSPAN1", "RNASE1", "CRYAB", "IFI27",
           "HLA-DQA1", "SCPEP1", "MUC15", "IL4I1", "ALDH1A")

salmon_quants_hg38.markers.emt <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% emt,]
rownames(salmon_quants_hg38.markers.emt) <- rowData(salmon_quants_hg38.markers.emt)$symbol
Heatmap(assay(salmon_quants_hg38.markers.emt, "izar"),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38$cluster.names,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE, col = getJetColors(circlize.cols = TRUE))

salmon_quants_hg38.markers.krt17 <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% krt17,]
rownames(salmon_quants_hg38.markers.krt17) <- rowData(salmon_quants_hg38.markers.krt17)$symbol
Heatmap(assay(salmon_quants_hg38.markers.krt17, "izar"),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38$cluster.names,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE, col = getJetColors(circlize.cols = TRUE))


## heatmap for Hui 2022-11-29
library(ComplexHeatmap)

secretory_markers <- read.delim("~/Downloads/secretory_markers.txt",
                                header=FALSE)
ciliated_markers <- read.delim("~/Downloads/ciliated_markers.txt",
                                header=FALSE)
## just look at secretory markers
secretory_markers <- c(secretory_markers$V1,
                       "MUC1")
## pull out just these markers
storm_all_pats.sec <- storm_all_pats[rowData(storm_all_pats)$symbol %in% secretory_markers,]

## pull out just secretory related clusters
storm_all_pats.sec <- storm_all_pats.sec[,grep("Secretory", storm_all_pats.sec$cluster.names)]

## make row names the gene name
rownames(storm_all_pats.sec) <- rowData(storm_all_pats.sec)$symbol

## order the secretory markers to be the same as Ian's heatmap
secretory_markers.order <- c("MUC1", "C3", "ASS1", "VTCN1",
                             "FAM107A", "MSLN", "KRT7", "ADAMTS1",
                             "PAX8", "EDN1", "LINC00937", "TTYH1",
                             "PLCB1", "FMOD", "PODXL", "ANO1",
                             "PKHD1L1", "OVGP1", "RSPO1")
secretory_markers.order.m <- match(secretory_markers.order,
                                   rownames(storm_all_pats.sec))
storm_all_pats.sec <- storm_all_pats.sec[secretory_markers.order.m,]

# careful
all(secretory_markers.order == rownames(storm_all_pats.sec))

assay(storm_all_pats.sec, "zscore") <- scale(assay(storm_all_pats.sec, "izar"))
# cap it
assay(storm_all_pats.sec, "zscore")[assay(storm_all_pats.sec, "zscore") > 3] <- 3
assay(storm_all_pats.sec, "zscore")[assay(storm_all_pats.sec, "zscore") < -3] <- -3

Heatmap(assay(storm_all_pats.sec, "zscore"),
        name = "Z-score",
        #column_split = storm_all_pats.sec$cluster.names,
        top_annotation = HeatmapAnnotation(Patient = storm_all_pats.sec$batch, col = list(Patient = c(patient1 = "lightgray", patient2 = "darkgray"))),
        show_column_names = FALSE,
        show_row_names = TRUE,
        #col = getJetColors(circlize.cols = TRUE),
        col = viridis(400, option = "D"),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        cluster_rows = FALSE)


