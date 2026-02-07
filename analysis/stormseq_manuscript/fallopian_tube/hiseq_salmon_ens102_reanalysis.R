## re-quant the hiseq patient 2 FTE
## load in the quants
## use velocessor
library(velocessor)

## load them in
quant.dirs <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/hiseq_patient2",
                         pattern = "quant", full.names = TRUE)
quants <- list.files(quant.dirs, pattern = "quant.sf", full.names = TRUE)
names(quants) <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/hiseq_patient2",
                            pattern = "quant")
salmon_quants_hg38.hiseq <- import_plate_txis(quants = quants,
                                        t2g = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.tx2gene.tsv",
                                        gtf = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.gtf")
saveRDS(salmon_quants_hg38.hiseq, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_raw_se.rds")

## izar transform to compare to the HGSOC data
salmon_quants_hg38.hiseq <- izar_transform(salmon_quants_hg38.hiseq)

## do PCA and usual cluster finding
#adjust for cell specific biases
salmon_quants_hg38.hiseq <- scran::computeSumFactors(salmon_quants_hg38.hiseq)
salmon_quants_hg38.hiseq <- scuttle::logNormCounts(salmon_quants_hg38.hiseq)

#fit a model for variance ~ expression
dec <- scran::modelGeneVar(salmon_quants_hg38.hiseq)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

#get the top 10% most variable genes
top.hvgs <- scran::getTopHVGs(dec, prop=0.1)
#PCA
salmon_quants_hg38.hiseq <- scater::runPCA(salmon_quants_hg38.hiseq, subset_row=top.hvgs,
                                     ncomponents = 20)

#find the number of PCs to retain
output <- scran::getClusteredPCs(reducedDim(salmon_quants_hg38.hiseq))
npcs <- metadata(output)$chosen
#15 PCs to hold onto
reducedDim(salmon_quants_hg38.hiseq, "PCAsub") <- reducedDim(salmon_quants_hg38.hiseq, "PCA")[,1:5,drop=FALSE]

#cluster using a shared nearest neighbor graph
g <- scran::buildSNNGraph(salmon_quants_hg38.hiseq, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership
salmon_quants_hg38.hiseq$cluster <- factor(cluster)
table(salmon_quants_hg38.hiseq$cluster)

#lib size
salmon_quants_hg38.hiseq$lib.size <- colSums(assay(salmon_quants_hg38.hiseq, "logcounts"))

## run density preserving t-SNE and UMAP
library(densvis)
set.seed(1988)

salmon_quants_hg38.hiseq.densne <- densne(reducedDim(salmon_quants_hg38.hiseq, "PCA"), dims = 2,
                                    verbose = TRUE, perplexity = 50)
salmon_quants_hg38.hiseq.densmap <- densmap(reducedDim(salmon_quants_hg38.hiseq, "PCA"), n_components = 2L,
                                      n_neighbors = 10L, metric = "euclidean")
salmon_quants_hg38.hiseq <- scater::runTSNE(salmon_quants_hg38.hiseq, dimred = "PCA")

## plot the den-SNE results
reducedDim(salmon_quants_hg38.hiseq, "denSNE") <- salmon_quants_hg38.hiseq.densne
reducedDim(salmon_quants_hg38.hiseq, "densMAP") <- salmon_quants_hg38.hiseq.densmap

## add some context
ens_genes <- get_ensembl_genes("102", "Homo sapiens")
ens_gene_match <- match(rownames(salmon_quants_hg38.hiseq), ens_genes$gene_id)
ens_genes <- ens_genes[ens_gene_match,]
rowData(salmon_quants_hg38.hiseq)$symbol <- ens_genes$symbol

colData(salmon_quants_hg38.hiseq)$epcam <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "EPCAM",]
colData(salmon_quants_hg38.hiseq)$cd31 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PECAM1",]
colData(salmon_quants_hg38.hiseq)$cd44 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD44",]
colData(salmon_quants_hg38.hiseq)$itga6 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ITGA6",]
#cluster 3 is likely our peg cells

colData(salmon_quants_hg38.hiseq)$cd34 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD34",]

colData(salmon_quants_hg38.hiseq)$tubb4 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "TUBB4A",]
colData(salmon_quants_hg38.hiseq)$pax8 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PAX8",]
colData(salmon_quants_hg38.hiseq)$cd45 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PTPRC",]
colData(salmon_quants_hg38.hiseq)$cd11c <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ITGAX",]
colData(salmon_quants_hg38.hiseq)$cd14 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD14",]
colData(salmon_quants_hg38.hiseq)$cd133 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PROM1",]

#usual stemmy markers
colData(salmon_quants_hg38.hiseq)$sox2 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "SOX2",]
colData(salmon_quants_hg38.hiseq)$nanog <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "NANOG",]
colData(salmon_quants_hg38.hiseq)$oct4 <- assays(salmon_quants_hg38.hiseq)$izar[rownames(salmon_quants_hg38.hiseq) == "ENSG00000204531",]

#Lan and Ron's paper on CA-MSCs
colData(salmon_quants_hg38.hiseq)$cd105 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ENG",]
colData(salmon_quants_hg38.hiseq)$cd90 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "THY1",]
colData(salmon_quants_hg38.hiseq)$cd73 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "NT5E",]

colData(salmon_quants_hg38.hiseq)$tgfb1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "TGFB1",]
colData(salmon_quants_hg38.hiseq)$bmp4 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "BMP4",]

#25 gene panel from Lan and Ron's stem cells paper
colData(salmon_quants_hg38.hiseq)$sox17 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "SOX17",]
colData(salmon_quants_hg38.hiseq)$crlf1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CRLF1",]

#ciliated markers from cancer cell 2020 paper
colData(salmon_quants_hg38.hiseq)$foxj1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "FOXJ1",]
colData(salmon_quants_hg38.hiseq)$ccdc17 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CCDC17",]
colData(salmon_quants_hg38.hiseq)$ccdc78 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CCDC78",]

#secretory markers from cancer cell 2020 paper
colData(salmon_quants_hg38.hiseq)$krt7 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "KRT7",]
#should be negative for CCDC17 and CD45

colData(salmon_quants_hg38.hiseq)$krt17 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "KRT17",]

#wt1
colData(salmon_quants_hg38.hiseq)$wt1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "WT1",]

colData(salmon_quants_hg38.hiseq)$cd3 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD3D",]
colData(salmon_quants_hg38.hiseq)$cd4 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD4",]
colData(salmon_quants_hg38.hiseq)$cd1a <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD1A",]
colData(salmon_quants_hg38.hiseq)$cd69 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CD69",]
colData(salmon_quants_hg38.hiseq)$cd103 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ITGAE",]
colData(salmon_quants_hg38.hiseq)$cd11b <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ITGAM",]
colData(salmon_quants_hg38.hiseq)$cd16b <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "FCGR3B",]

colData(salmon_quants_hg38.hiseq)$calb2 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "CALB2",]
colData(salmon_quants_hg38.hiseq)$tp53 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "TP53",]
colData(salmon_quants_hg38.hiseq)$klf4 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "KLF4",]

colData(salmon_quants_hg38.hiseq)$sparcl1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "SPARCL1",]
colData(salmon_quants_hg38.hiseq)$rbp1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "RBP1",]
colData(salmon_quants_hg38.hiseq)$aldh1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ALDH1A1",]

colData(salmon_quants_hg38.hiseq)$mcam <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "MCAM",]
colData(salmon_quants_hg38.hiseq)$pdgfrb <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PDGFRB",]
colData(salmon_quants_hg38.hiseq)$susd2 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "SUSD2",]
colData(salmon_quants_hg38.hiseq)$runx3 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "RUNX3",]
colData(salmon_quants_hg38.hiseq)$igfbp5 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "IGFBP5",]
colData(salmon_quants_hg38.hiseq)$acta2 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "ACTA2",]
colData(salmon_quants_hg38.hiseq)$igfbp3 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "IGFBP3",]
colData(salmon_quants_hg38.hiseq)$pecam1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "PECAM1",]
colData(salmon_quants_hg38.hiseq)$ovgp1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "OVGP1",]
colData(salmon_quants_hg38.hiseq)$jam3 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "JAM3",]
colData(salmon_quants_hg38.hiseq)$sparcl1 <- assays(salmon_quants_hg38.hiseq)$izar[rowData(salmon_quants_hg38.hiseq)$symbol == "SPARCL1",]


## looks great - use the densMAP embedding
scater::plotTSNE(salmon_quants_hg38.hiseq, colour_by = "ovgp1")
scater::plotReducedDim(salmon_quants_hg38.hiseq, "denSNE", colour_by = "cluster")
scater::plotReducedDim(salmon_quants_hg38.hiseq, "densMAP", colour_by = "epcam", by_exprs_values = "izar")

saveRDS(salmon_quants_hg38.hiseq, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")

## strip down for import into Python
salmon_quants_hg38.stripped <- salmon_quants_hg38.hiseq
metadata(salmon_quants_hg38.stripped) <- list() ## zero out metadata
assays(salmon_quants_hg38.stripped) <- list(
  counts = round(assay(salmon_quants_hg38.hiseq, "spliced")),
  spliced = round(assay(salmon_quants_hg38.hiseq, "spliced")),
  unspliced = round(assay(salmon_quants_hg38.hiseq, "unspliced"))
)

rownames(salmon_quants_hg38.stripped) <- rowData(salmon_quants_hg38.stripped)$symbol
rownames(salmon_quants_hg38.stripped) <- make.unique(rownames(salmon_quants_hg38.stripped), sep = "_")
saveRDS(salmon_quants_hg38.stripped, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings_forscvelo.rds")
library(zellkonverter)
writeH5AD(salmon_quants_hg38.stripped, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings_forscvelo.h5ad")

#roll my own plotting function...
library(ggplot2)
library(cowplot)

to_plot <- as.data.frame(reducedDim(salmon_quants_hg38.hiseq, "TSNE"))
colnames(to_plot) <- c("Dim_1", "Dim_2")

epcam <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38.hiseq)$epcam)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EpCAM gene expression") +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

cd31 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38.hiseq)$cd31)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("CD31 gene expression") +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

#cd105 cd90 cd73
stem_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38.hiseq)[,c("cd105", "cd90", "cd73")]))
msc <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = stem_sum_logcounts)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("MSC-like marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

#secretory cells
secretory_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38.hiseq)[,c("pax8", "krt7")]))
secretory <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = secretory_sum_logcounts)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
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

#ciliated
ciliated_sum_logcounts <- rowMeans(as.data.frame(colData(salmon_quants_hg38.hiseq)[,c("foxj1", "ccdc17", "ccdc78")]))
ciliated <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = ciliated_sum_logcounts)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
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

#monocyte
monocyte_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38.hiseq)[,c("cd45", "cd11c", "cd14")]))
monocyte <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = monocyte_sum_logcounts)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
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

#CD of interest
cd_ofinterest <- rowSums(as.data.frame(colData(salmon_quants_hg38.hiseq)[,c("cd34")]))
cd_ofinterest_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = cd_ofinterest)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("CD34 marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(epcam, msc, secretory, ciliated, monocyte, cd_ofinterest_plot, labels = "AUTO", nrow = 3, ncol = 2)

plot_grid(epcam, secretory, ciliated, monocyte, nrow = 2, ncol = 2)

## look at marker gene expression?
load('~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/novaseq_marker_genes_patient1.rda')

salmon_quants_hg38.hiseq.markers <- salmon_quants_hg38.hiseq[rowData(salmon_quants_hg38.hiseq)$symbol %in% marker_genes_by_cluster$marker_genes,]

library(dplyr)
salmon_quants_hg38$cluster.collapse <- salmon_quants_hg38$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "Ciliated",
    . == 2 ~ "Secretory",
    . == 3 ~ "MSC_like",
    . == 4 ~ "Ciliated",
    . == 5 ~ "Branch",
    . == 6 ~ "Immune",
    . == 7 ~ "Secretory",
    . == 8 ~ "Secretory"
  ))
salmon_quants_hg38$cluster.collapse <- salmon_quants_hg38$cluster.collapse$cluster_cat



## add in the surface level expression
#let's tack on the cell surface expression of EpCAM from flow
#thanks Rachael!
surface_expression <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/figures/384_well_sort_SI_Events_converted_hiseq.txt",
                                 stringsAsFactors = FALSE)
rownames(surface_expression) <- paste0(surface_expression$Sort.Index.Y,
                                       surface_expression$Sort.Index.X)
epcam_surface_expression <- data.frame(epcam = surface_expression$EpCAM,
                                       ssc_a = surface_expression$SSC.A,
                                       ssc_w = surface_expression$SSC.W,
                                       ssc_h = surface_expression$SSC.H)
rownames(epcam_surface_expression) <- rownames(surface_expression)
epcam_surface_expression$well_id <- rownames(epcam_surface_expression)
#match em up
# first need to bring in and parse the well IDs.... argh.
barcode_assoc_1 <- read.delim("~/Downloads/barcodeAssociationTable_1.txt",
                              sep = ",")
barcode_assoc_2 <- read.delim("~/Downloads/barcodeAssociationTable_2.txt",
                              sep = ",")
barcode_assoc_1$short_name <- gsub(".*-", "", barcode_assoc_1$AccessionID)
barcode_assoc_2$short_name <- gsub(".*-", "", barcode_assoc_2$AccessionID)
all_barcodes <- rbind(barcode_assoc_1[,c(8,2)],
                      barcode_assoc_2[,c(8,2)])

# helper
rmAfterNthDelim <- function(vec, delim="-", n=3) {
  return(substr(vec, 1, sapply(gregexpr(delim, vec), "[", n) - 1))
}

salmon_quants_hg38.hiseq$short_name <- rmAfterNthDelim(salmon_quants_hg38.hiseq$sample,
                                                       delim="_" ,n = 2)
salmon_quants_hg38.hiseq$short_name <- gsub(".*_", "", salmon_quants_hg38.hiseq$short_name)

## restrict the epcam surface expressin matrix to only those found in the counts obj
short_name_match <- match(salmon_quants_hg38.hiseq$short_name,
                          all_barcodes$short_name)
all_barcodes <- all_barcodes[short_name_match,]
#all(all_barcodes$short_name == salmon_quants_hg38.hiseq$short_name)
#TRUE
salmon_quants_hg38.hiseq$clean.sample.id <- all_barcodes$ClientAccessionID
salmon_quants_hg38.hiseq <- salmon_quants_hg38.hiseq[,salmon_quants_hg38.hiseq$clean.sample.id %in% rownames(epcam_surface_expression)]
## clean up sample names to wells
ese.match <- match(rownames(epcam_surface_expression), salmon_quants_hg38.hiseq$clean.sample.id)
ese.match.clean <- ese.match[!is.na(ese.match)]
epcam_surface_expression$match <- ese.match
epcam_surface_expression <- epcam_surface_expression[!is.na(epcam_surface_expression$match),]
salmon_quants_hg38.hiseq <- salmon_quants_hg38.hiseq[,ese.match.clean]
all(epcam_surface_expression$well_id == salmon_quants_hg38.hiseq$clean.sample.id)
#TRUE
#epcam_surface_expression <- subset(epcam_surface_expression, !is.na(epcam_surface_expression$epcam))
colData(salmon_quants_hg38.hiseq)$epcam_surface_expression <- epcam_surface_expression$epcam
colData(salmon_quants_hg38.hiseq)$ssc_a <- epcam_surface_expression$ssc_a
colData(salmon_quants_hg38.hiseq)$ssc_h <- epcam_surface_expression$ssc_h
colData(salmon_quants_hg38.hiseq)$ssc_w <- epcam_surface_expression$ssc_w

quantiles_epcam <- gtools::quantcut(salmon_quants_hg38.hiseq$epcam_surface_expression)
quantiles_epcam <- forcats::fct_recode(quantiles_epcam,
                                       `0-25%` = "[616,2.99e+03]",
                                       `25-50%` = "(2.99e+03,4.45e+03]",
                                       `50-75%` = "(4.45e+03,5.88e+03]",
                                       `75-100%` = "(5.88e+03,2.86e+04]")
salmon_quants_hg38.hiseq$epcam_surface_expression_quantiles <- quantiles_epcam

## now let's call marker genes for clusters
scater::plotReducedDim(salmon_quants_hg38.hiseq, "TSNE", colour_by = "cluster")
scater::plotReducedDim(salmon_quants_hg38.hiseq, "TSNE", colour_by = "epcam_surface_expression")
scater::plotReducedDim(salmon_quants_hg38.hiseq, "TSNE", colour_by = "epcam_surface_expression_quantiles")

## some ggplots to look at putative doublets and color by cluster
to_plot <- as.data.frame(colData(salmon_quants_hg38.hiseq))
#to_plot_sub <- to_plot[to_plot$cluster == 3,]
ssc_plot_a <- ggplot(to_plot, aes(x = ssc_a, y = ssc_h,
                                  color = cluster)) +
  geom_point() +
  scale_color_viridis_d(begin = 0.2, name = "Cluster") +
  # geom_point(data = to_plot_sub,
  #            aes(x = ssc_a, y = ssc_h, color = "red")) +
  #scale_color_manual(values = to_plot$color) +
  theme_half_open() +
  ylab("SSC-H") +
  xlab("SSC-A") +
  ylim(c(0, 140)) +
  xlim(c(0, 140)) +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  #ggtitle("SSC-A vs SSC-H colored by cluster") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "None",
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ssc_plot_b <- ggplot(to_plot, aes(x = ssc_h, y = ssc_w,
                                  color = cluster)) +
  geom_point() +
  #geom_point(alpha = 0.5) +
  # geom_point(data = to_plot_sub,
  #            aes(x = ssc_h, y = ssc_w, color = "red")) +
  scale_color_viridis_d(begin = 0.2, name = "Cluster") +
  theme_half_open() +
  ylab("SSC-W") +
  xlab("SSC-H") +
  ylim(c(0, 140)) +
  xlim(c(0, 140)) +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  #ggtitle("SSC-H vs SSC-W colored by cluster") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    #legend.position = "None",
    legend.text = element_text(size = 16)
    #legend.title = element_blank()
  )

plot_grid(ssc_plot_a, ssc_plot_b, ncol = 2, nrow = 1, rel_widths = c(0.8,1))

## plot the quartiles of epcam expression
to_plot$Dim1 <- reducedDim(salmon_quants_hg38.hiseq, "TSNE")[,1]
to_plot$Dim2 <- reducedDim(salmon_quants_hg38.hiseq, "TSNE")[,2]
epcam_quantiles_plot <- ggplot(to_plot, aes(x = Dim1, y = Dim2,
                                  color = epcam_surface_expression_quantiles)) +
  geom_point() +
  scale_color_brewer(palette = "Set1", name = "EpCAM surface\nexpression quantiles") +
  theme_half_open() +
  ylab("Dimension 1") +
  xlab("Dimension 2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(hjust = 0.5, size = 16)
  )

epcam_logexp_plot <- ggplot(to_plot, aes(x = Dim1, y = Dim2,
                                            color = log2(epcam_surface_expression))) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "EpCAM log2\nsurface expression") +
  theme_half_open() +
  ylab("Dimension 1") +
  xlab("Dimension 2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(hjust = 0.5, size = 16)
  )

hiseq_clust_plot <- ggplot(to_plot, aes(x = Dim1, y = Dim2,
                                         color = cluster)) +
  geom_point() +
  scale_color_viridis_d(begin = 0.2, name = "Cluster") +
  theme_half_open() +
  ylab("Dimension 1") +
  xlab("Dimension 2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    #legend.position = "None",
    legend.text = element_text(size = 16)
    #legend.title = element_blank()
  )
