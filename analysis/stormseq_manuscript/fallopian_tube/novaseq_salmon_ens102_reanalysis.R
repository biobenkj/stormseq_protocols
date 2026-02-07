## Re-quantify the 20M sub-sampled files to GRCh38 release 102
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
})

## following something similar to: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
## modifications after talking with Charlotte, Michael, Avi, and Rob

gtf <- "Homo_sapiens.GRCh38.102.gtf.gz"

## We want to extract the spliced and unspliced tx instead of 
## spliced and introns since we have full gene structure with smart-seq

grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "unspliced"),
  verbose = TRUE
)
# Import genomic features from the file as a GRanges object ... OK
# Prepare the 'metadata' data frame ... OK
# Make the TxDb object ... OK
# 'select()' returned 1:1 mapping between keys and columns
# Extracting spliced transcript features
# Extracting unspliced transcript features
# Warning message:
# In .get_cds_IDX(mcols0$type, mcols0$phase) :
# The "phase" metadata column contains non-NA values for features of type
# stop_codon. This information was ignored.

## On to extracting things from the reference genome fasta
## this isn't particularly happy with a symlink on mounted storage
## specify a relative path to the actual file if using mounted storage
genome <- Biostrings::readDNAStringSet("Homo_sapiens.GRCh38.dna.primary_assembly.fa")
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)
Biostrings::writeXStringSet(
  seqs, filepath = "Homo_sapiens.GRCh38.dna.primary_assembly.expanded.fa"
)

## export the annotations
eisaR::exportToGtf(
  grl, 
  filepath = "Homo_sapiens.GRCh38.102.expanded.gtf"
)

## write out the feature names for the index
## namely the spliced and unspliced transcripts for index generation

head(metadata(grl)$corrgene)
# spliced         unspliced
# 1 ENSG00000223972 ENSG00000223972-U
# 2 ENSG00000243485 ENSG00000243485-U
# 3 ENSG00000284332 ENSG00000284332-U
# 4 ENSG00000268020 ENSG00000268020-U
# 5 ENSG00000240361 ENSG00000240361-U
# 6 ENSG00000186092 ENSG00000186092-U

write.table(
  metadata(grl)$corrgene, 
  file = "Homo_sapiens.GRCh38.102.expanded.features.tsv",
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

## transcript (spliced and unspliced) map to genes
df <- eisaR::getTx2Gene(
  grl, filepath = "Homo_sapiens.GRCh38.102.expanded.tx2gene.tsv"
)

## now make the linked transcriptome map with tximeta
tximeta::makeLinkedTxome(
  indexDir = "Homo_sapiens.GRCh38.102.genome.expanded.sidx", 
  source = "Ensembl", genome = "GRCh38", 
  organism = "Homo sapiens", release = "102", 
  fasta = "Homo_sapiens.GRCh38.dna.primary_assembly.expanded.fa", 
  gtf = "Homo_sapiens.GRCh38.102.expanded.gtf", 
  write = TRUE, jsonFile = "Homo_sapiens.GRCh38.102.expanded.json"
)

## load in the quants
## use velocessor
library(velocessor)

## load them in
quant.dirs <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant",
                         pattern = "_quant", full.names = TRUE)
quants <- list.files(quant.dirs, pattern = "quant.sf", full.names = TRUE)
names(quants) <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant",
                            pattern = "_quant")
salmon_quants_hg38 <- import_plate_txis(quants = quants,
                                        t2g = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.tx2gene.tsv",
                                        gtf = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.gtf")
saveRDS(salmon_quants_hg38, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_raw_se.rds")

## izar transform to compare to the HGSOC data
salmon_quants_hg38 <- izar_transform(salmon_quants_hg38)

## do PCA and usual cluster finding
#adjust for cell specific biases
salmon_quants_hg38 <- scran::computeSumFactors(salmon_quants_hg38)
salmon_quants_hg38 <- scuttle::logNormCounts(salmon_quants_hg38)

#fit a model for variance ~ expression
dec <- scran::modelGeneVar(salmon_quants_hg38)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

#get the top 10% most variable genes
top.hvgs <- scran::getTopHVGs(dec, prop=0.1)
#PCA
salmon_quants_hg38 <- scater::runPCA(salmon_quants_hg38, subset_row=top.hvgs,
                                     ncomponents = 20)

#find the number of PCs to retain
output <- scran::getClusteredPCs(reducedDim(salmon_quants_hg38))
npcs <- metadata(output)$chosen
#7 PCs to hold onto
reducedDim(salmon_quants_hg38, "PCAsub") <- reducedDim(salmon_quants_hg38, "PCA")[,1:npcs,drop=FALSE]

#cluster using a shared nearest neighbor graph
g <- scran::buildSNNGraph(salmon_quants_hg38, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership
salmon_quants_hg38$cluster <- factor(cluster)
table(salmon_quants_hg38$cluster)

#lib size
salmon_quants_hg38$lib.size <- colSums(assay(salmon_quants_hg38, "logcounts"))

## run density preserving t-SNE and UMAP
library(densvis)
set.seed(1988)

salmon_quants_hg38.densne <- densne(reducedDim(salmon_quants_hg38, "PCA"), dims = 3,
                                    verbose = TRUE, perplexity = 50)
salmon_quants_hg38.densmap <- densmap(reducedDim(salmon_quants_hg38, "PCA"), n_components = 3L,
                                      n_neighbors = 15L, metric = "euclidean")
salmon_quants_hg38 <- scater::runTSNE(salmon_quants_hg38, dimred = "PCA")

## plot the den-SNE results
reducedDim(salmon_quants_hg38, "denSNE") <- salmon_quants_hg38.densne
reducedDim(salmon_quants_hg38, "densMAP") <- salmon_quants_hg38.densmap

## looks great - use the densMAP embedding
scater::plotTSNE(salmon_quants_hg38, colour_by = "cluster")
scater::plotReducedDim(salmon_quants_hg38, "denSNE", colour_by = "cluster")
scater::plotReducedDim(salmon_quants_hg38, "densMAP", colour_by = "cluster", ncomponents = c(1,2))
scater::plotReducedDim(salmon_quants_hg38, "densMAP", colour_by = "cluster", ncomponents = c(2,3))

saveRDS(salmon_quants_hg38, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")

## strip down for import into Python
salmon_quants_hg38.stripped <- salmon_quants_hg38
metadata(salmon_quants_hg38.stripped) <- list() ## zero out metadata
assays(salmon_quants_hg38.stripped) <- list(
  counts = round(assay(salmon_quants_hg38, "spliced")),
  spliced = round(assay(salmon_quants_hg38, "spliced")),
  unspliced = round(assay(salmon_quants_hg38, "unspliced"))
)

salmon_quants_hg38.stripped <- salmon_quants_hg38.stripped[,!salmon_quants_hg38.stripped$cluster %in% 6]
salmon_quants_hg38.stripped$cluster <- factor(salmon_quants_hg38.stripped$cluster, levels = c(1,2,3,4,5,7,8))
rownames(salmon_quants_hg38.stripped) <- rowData(salmon_quants_hg38.stripped)$symbol
rownames(salmon_quants_hg38.stripped) <- make.unique(rownames(salmon_quants_hg38.stripped), sep = "_")
saveRDS(salmon_quants_hg38.stripped, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings_forscvelo.rds")
library(zellkonverter)
writeH5AD(salmon_quants_hg38.stripped, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings_forscvelo.h5ad")
## perhaps worth plotting in PCA space as well
scater::plotPCA(salmon_quants_hg38, colour_by = "cluster")
scater::plotPCA(salmon_quants_hg38, colour_by = "cluster", ncomponents = c(1,2))
## whoa! in PCA space, we get a similar trajectory!
## the perk of PCA space is we can _directly_ examine the loadings

## pull off the top genes driving PC1
fte_pca <- attr(reducedDim(salmon_quants_hg38, "PCA"), "rotation")
fte_pca.pc1 <- fte_pca[order(abs(fte_pca[,1]), decreasing = TRUE),]
fte_pca.pc1.top100 <- rownames(fte_pca.pc1)[1:100]
write.table(fte_pca.pc1.top100, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/top100_ensgenes_pc1_nova_fte_downsampled_grch38.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# fte_pca.pc1.top200 <- rownames(fte_pca.pc1)[1:200]
# write.table(fte_pca.pc1.top200, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/top200_ensgenes_pc1_nova_fte_downsampled_grch38.txt",
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## neat! the top genes are for cilia formation and movement.
## there's also mention of nuclear lamen changes - chromatin remodeling

## make a quick bar chart from results


## color some marker genes
#color by marker gene expression
colData(salmon_quants_hg38)$epcam <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "EPCAM",]
colData(salmon_quants_hg38)$cd44 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD44",]
colData(salmon_quants_hg38)$itga6 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ITGA6",]
#cluster 3 is likely our peg cells

colData(salmon_quants_hg38)$cd34 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD34",]

colData(salmon_quants_hg38)$tubb4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "TUBB4B",]
colData(salmon_quants_hg38)$pax8 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PAX8",]
colData(salmon_quants_hg38)$cd45 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PTPRC",]
colData(salmon_quants_hg38)$cd11c <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ITGAX",]
colData(salmon_quants_hg38)$cd14 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD14",]
colData(salmon_quants_hg38)$cd133 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PROM1",]

#usual stemmy markers
colData(salmon_quants_hg38)$sox2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SOX2",]
colData(salmon_quants_hg38)$nanog <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "NANOG",]
colData(salmon_quants_hg38)$oct4 <- assays(salmon_quants_hg38)$izar[rownames(salmon_quants_hg38) == "ENSG00000204531",]

#Lan and Ron's paper on CA-MSCs
colData(salmon_quants_hg38)$cd105 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ENG",]
colData(salmon_quants_hg38)$cd90 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "THY1",]
colData(salmon_quants_hg38)$cd73 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "NT5E",]

colData(salmon_quants_hg38)$tgfb1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "TGFB1",]
colData(salmon_quants_hg38)$bmp4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "BMP4",]

#25 gene panel from Lan and Ron's stem cells paper
colData(salmon_quants_hg38)$sox17 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SOX17",]
colData(salmon_quants_hg38)$crlf1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CRLF1",]

#ciliated markers from cancer cell 2020 paper
colData(salmon_quants_hg38)$foxj1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FOXJ1",]
colData(salmon_quants_hg38)$ccdc17 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CCDC17",]
colData(salmon_quants_hg38)$ccdc78 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CCDC78",]

#secretory markers from cancer cell 2020 paper
colData(salmon_quants_hg38)$krt7 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "KRT7",]
#should be negative for CCDC17 and CD45

colData(salmon_quants_hg38)$krt17 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "KRT17",]

#wt1
colData(salmon_quants_hg38)$wt1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "WT1",]

colData(salmon_quants_hg38)$cd3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD3D",]
colData(salmon_quants_hg38)$cd4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD4",]
colData(salmon_quants_hg38)$cd1a <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD1A",]
colData(salmon_quants_hg38)$cd69 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD69",]
colData(salmon_quants_hg38)$cd103 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ITGAE",]
colData(salmon_quants_hg38)$cd11b <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ITGAM",]
colData(salmon_quants_hg38)$cd16b <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FCGR3B",]

colData(salmon_quants_hg38)$calb2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CALB2",]
colData(salmon_quants_hg38)$tp53 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "TP53",]
colData(salmon_quants_hg38)$klf4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "KLF4",]

colData(salmon_quants_hg38)$sparcl1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SPARCL1",]
colData(salmon_quants_hg38)$rbp1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "RBP1",]
colData(salmon_quants_hg38)$aldh1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ALDH1A1",]
colData(salmon_quants_hg38)$aldh1a3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ALDH1A3",]

colData(salmon_quants_hg38)$mcam <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MCAM",]
colData(salmon_quants_hg38)$pdgfrb <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PDGFRB",]
colData(salmon_quants_hg38)$susd2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SUSD2",]
colData(salmon_quants_hg38)$runx3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "RUNX3",]
colData(salmon_quants_hg38)$igfbp5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IGFBP5",]
colData(salmon_quants_hg38)$acta2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ACTA2",]
colData(salmon_quants_hg38)$igfbp3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IGFBP3",]
colData(salmon_quants_hg38)$pecam1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PECAM1",]
colData(salmon_quants_hg38)$ovgp1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "OVGP1",]
colData(salmon_quants_hg38)$jam3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "JAM3",]
colData(salmon_quants_hg38)$sparcl1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SPARCL1",]

colData(salmon_quants_hg38)$cldn5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CLDN5",]
colData(salmon_quants_hg38)$fabp4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FABP4",]
colData(salmon_quants_hg38)$pprc1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PPRC1",]

## stromal markers from hu et al. 2020
colData(salmon_quants_hg38)$col1a2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "COL1A2",]
colData(salmon_quants_hg38)$col3a1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "COL3A1",]

colData(salmon_quants_hg38)$msln <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MSLN",]
colData(salmon_quants_hg38)$esr1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ESR1",]
colData(salmon_quants_hg38)$brca1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "BRCA1",]

## find them endothelial cells
colData(salmon_quants_hg38)$vwf <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "VWF",]
colData(salmon_quants_hg38)$cdh5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CDH5",]
colData(salmon_quants_hg38)$flt1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FLT1",]
colData(salmon_quants_hg38)$kdr <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "KDR",]
colData(salmon_quants_hg38)$pecam1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PECAM1",]
colData(salmon_quants_hg38)$cd34 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "CD34",]

colData(salmon_quants_hg38)$tmem173 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "STING1",]

colData(salmon_quants_hg38)$sema4a <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SEMA4A",]


colData(salmon_quants_hg38)$gata2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "GATA2",]

colData(salmon_quants_hg38)$angpt2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ANGPT2",]

colData(salmon_quants_hg38)$mir181a1hg <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MIR181A1HG",]
colData(salmon_quants_hg38)$mir200chg <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MIR200CHG",]

## compare to Jun's data
colData(salmon_quants_hg38)$lgr5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "LGR5",]
colData(salmon_quants_hg38)$prrx1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PRRX1",]
colData(salmon_quants_hg38)$pgr <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PGR",]

#roll my own plotting function...
library(ggplot2)
library(cowplot)

to_plot <- as.data.frame(reducedDim(salmon_quants_hg38, "densMAP"))
colnames(to_plot) <- c("Dim_1", "Dim_2")

epcam_cs <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = log2(salmon_quants_hg38$epcam_surface_expression))) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EpCAM cell surface expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

epcam <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$epcam)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
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
cd31 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$pecam1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("CD31 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )
plot_grid(epcam_cs, epcam, cd31,
          ncol = 3, nrow = 1, align = "hv")

clst <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = salmon_quants_hg38$cluster)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_viridis_d() +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  #ggtitle("FTE Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "None"
  )
clst
endo_tpm <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd34", "vwf", "cdh5", "pecam1",
                                                                 "kdr", "flt1")]))
endothelial <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = endo_tpm)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Endothelial cell gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

jun_epi_pro_pop1 <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("lgr5", "pgr")]))
jun_epi_pro_pop1_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = jun_epi_pro_pop1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("LGR5+/PGR+ gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

jun_epi_pro_pop2 <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("acta2", "prrx1")]))
jun_epi_pro_pop2_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = jun_epi_pro_pop2)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("ACTA2+/PRRX1+ gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(jun_epi_pro_pop1_plot, jun_epi_pro_pop2_plot, labels = "AUTO",
          ncol = 2, nrow = 1, align = "hv")

lgr5 <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("lgr5")]))
lgr5_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = lgr5)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("LGR5 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

pgr <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("pgr")]))
pgr_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = pgr)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("PGR gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

acta2 <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("acta2")]))
acta2_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = acta2)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("ACTA2 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

prrx1 <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("prrx1")]))
prrx1_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = prrx1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("PRRX1 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(lgr5_plot, pgr_plot, acta2_plot, prrx1_plot, labels = "AUTO",
          ncol = 2, nrow = 2, align = "hv")

to_plot <- as.data.frame(reducedDim(salmon_quants_hg38, "densMAP"))
colnames(to_plot) <- c("Dim_1", "Dim_2")
telocyte_tpm <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd117", "cd34", "des", "gja1", "pdgfrb")]))
telocyte <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = telocyte_tpm)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Telocyte gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

#abcg2
#slc2a3
#cd34

#cd105
#cd34
#slc2a3


tz <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("abcg2", "slc2a3", "cd34")]))
tz_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = tz)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("TZ root-like marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

bkj <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd105", "slc2a3", "cd34")]))
bkj_plot <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = bkj)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("BKJ root-like marker gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

#cd105 cd90 cd73
stem_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd105", "cd90", "cd73")]))
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
secretory_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("pax8", "krt7")]))
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
ciliated_sum_logcounts <- rowMeans(as.data.frame(colData(salmon_quants_hg38)[,c("foxj1", "ccdc17", "ccdc78")]))
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
monocyte_sum_logcounts <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd45", "cd11c", "cd14")]))
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
cd_ofinterest <- rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("cd34")]))
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

## goooood

## look at interferons
colData(salmon_quants_hg38)$ifnb1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNB1",]

colData(salmon_quants_hg38)$ifna1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA1",]
colData(salmon_quants_hg38)$ifna2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA2",]
colData(salmon_quants_hg38)$ifna4 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA4",]
colData(salmon_quants_hg38)$ifna5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA5",]
colData(salmon_quants_hg38)$ifna6 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA6",]
colData(salmon_quants_hg38)$ifna7 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNA7",]

# epithelial co-regulated interferons?
colData(salmon_quants_hg38)$ifnl1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNL1",]
colData(salmon_quants_hg38)$ifnl2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNL2",]
colData(salmon_quants_hg38)$ifnl3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFNL3",]

# interferon regulation
colData(salmon_quants_hg38)$irf3 <- assays(salmon_quants_hg38)$logcounts[rowData(salmon_quants_hg38)$symbol == "IRF3",]
colData(salmon_quants_hg38)$irf7 <- assays(salmon_quants_hg38)$logcounts[rowData(salmon_quants_hg38)$symbol == "IRF7",]
# more IRF3 expression than IRF7

# RNA editing
colData(salmon_quants_hg38)$adar1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ADAR",]
# definitely ADAR1 expression in both arms, but higher (very slightly) in ciliated

# ISG regulators
colData(salmon_quants_hg38)$mx1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MX1",]
colData(salmon_quants_hg38)$mx2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MX2",]

colData(salmon_quants_hg38)$ifit1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "IFIT1",]
colData(salmon_quants_hg38)$oas1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "OAS1",]

# downstream JAK/STAT signaling
colData(salmon_quants_hg38)$stat1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "STAT1",]
colData(salmon_quants_hg38)$tbk1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "TBK1",]
# zinc fingers
colData(salmon_quants_hg38)$znf382 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ZNF382",]
colData(salmon_quants_hg38)$znf384 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ZNF384",]

# testing mouse model hypothesis
colData(salmon_quants_hg38)$krt5 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "KRT5",]
colData(salmon_quants_hg38)$tp73 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "TP73",]
colData(salmon_quants_hg38)$cd133 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "PROM1",]
colData(salmon_quants_hg38)$slc1a3 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "SLC1A3",]

ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("krt5", "tp73", "cd133")])))) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("KRT5/TP73/CD133 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("tp53", "slc1a3", "pax8")])))) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("SLC1A3/TP53/PAX8 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = rowSums(as.data.frame(colData(salmon_quants_hg38)[,c("pax8")])))) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("PAX8 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

## add in the surface level expression
#let's tack on the cell surface expression of EpCAM from flow
#thanks Rachael!
surface_expression <- read.table("~/Documents/manuscripts/metabolic_flow_fte_2020/figures/Converted384wellIndexSort.txt",
                                 stringsAsFactors = F, header = T)
rownames(surface_expression) <- paste0(surface_expression$Row,
                                       surface_expression$Column)
epcam_surface_expression <- data.frame(epcam = surface_expression$EpCAM,
                                       ssc_a = surface_expression$SSC.A,
                                       ssc_w = surface_expression$SSC.W,
                                       ssc_h = surface_expression$SSC.H)
rownames(epcam_surface_expression) <- rownames(surface_expression)
epcam_surface_expression$well_id <- rownames(epcam_surface_expression)
#match em up
## clean up sample names to wells
salmon_quants_hg38$clean.sample.id <- gsub("_quant", "", salmon_quants_hg38$sample)
salmon_quants_hg38$clean.sample.id <- gsub("\\_.*$", "", salmon_quants_hg38$clean.sample.id)
ese.match <- match(salmon_quants_hg38$clean.sample.id, rownames(epcam_surface_expression))
epcam_surface_expression <- as.data.frame(epcam_surface_expression[ese.match,])
#epcam_surface_expression <- subset(epcam_surface_expression, !is.na(epcam_surface_expression$epcam))
colData(salmon_quants_hg38)$epcam_surface_expression <- epcam_surface_expression$epcam
colData(salmon_quants_hg38)$ssc_a <- epcam_surface_expression$ssc_a
colData(salmon_quants_hg38)$ssc_h <- epcam_surface_expression$ssc_h
colData(salmon_quants_hg38)$ssc_w <- epcam_surface_expression$ssc_w

quantiles_epcam <- gtools::quantcut(salmon_quants_hg38$epcam_surface_expression)
quantiles_epcam <- forcats::fct_recode(quantiles_epcam,
                                       `0-25%` = "[479,1.43e+03]",
                                       `25-50%` = "(1.43e+03,2.2e+03]",
                                       `50-75%` = "(2.2e+03,3.14e+03]",
                                       `75-100%` = "(3.14e+03,1.98e+04]")
salmon_quants_hg38$epcam_surface_expression_quantiles <- quantiles_epcam

## now let's call marker genes for clusters
scater::plotReducedDim(salmon_quants_hg38, "densMAP", colour_by = "cluster")
scater::plotReducedDim(salmon_quants_hg38, "densMAP", colour_by = "epcam_surface_expression")
scater::plotReducedDim(salmon_quants_hg38, "densMAP", colour_by = "epcam_surface_expression_quantiles")

## some ggplots to look at putative doublets and color by cluster
to_plot <- as.data.frame(colData(salmon_quants_hg38))
to_plot_sub <- to_plot[to_plot$cluster == 3,]
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
    legend.text = element_text(size = 16),
    #legend.title = element_blank()
  )

plot_grid(ssc_plot_a, ssc_plot_b, ncol = 2, nrow = 1, rel_widths = c(0.8,1))

## plot the quartiles of epcam expression
to_plot$Dim1 <- reducedDim(salmon_quants_hg38, "TSNE")[,1]
to_plot$Dim2 <- reducedDim(salmon_quants_hg38, "TSNE")[,2]
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

novaseq_clust_plot <- ggplot(to_plot, aes(x = Dim1, y = Dim2,
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

to_plot <- as.data.frame(reducedDim(salmon_quants_hg38, "densMAP"))
colnames(to_plot) <- c("Dim_1", "Dim_2")

epcam <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$epcam)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
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
epcam_surface <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = salmon_quants_hg38$epcam_surface_expression)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("EpCAM surface expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(epcam, epcam_surface, labels = "AUTO", nrow = 1, ncol = 2)


## mask off immune cells
salmon_quants_hg38.mask <- salmon_quants_hg38[,!salmon_quants_hg38$cluster %in% 6]
salmon_quants_hg38.mask$cluster <- factor(salmon_quants_hg38$cluster, levels = c(1,2,3,4,5,7,8))

## collapse clusters
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

salmon_quants_hg38.sub <- salmon_quants_hg38[,salmon_quants_hg38$cluster.collapse %in% c("Ciliated",
                                                                                         "Secretory")]

markers.sce.cluster <- scran::findMarkers(salmon_quants_hg38, salmon_quants_hg38$cluster.collapse,
                                   test = "wilcox", direction = "up",
                                   #gene.names = rowData(salmon_quants_hg38)$symbol,
                                   pval.type = "all", assay.type = "izar")

markers.sce.cluster.nocollapse <- scran::findMarkers(salmon_quants_hg38, salmon_quants_hg38$cluster,
                                          test = "wilcox", direction = "up",
                                          gene.names = rowData(salmon_quants_hg38)$symbol,
                                          pval.type = "some", assay.type = "izar")
#cluster_annots <- data.frame(AUC_clust = factor(c("ciliated_arm", "secretory_arm", "msc_cluster",
#                                                  "ciliated_cluster", "arm_center", "immune_cluster",
#                                                  "secretory_clusterA", "secretory_clusterB")))
#rownames(cluster_annots) <- paste0("AUC.", 1:8)
#cluster_cols <- list(AUC_clust = RColorBrewer::brewer.pal(8, "PuBu"))
#names(cluster_cols$AUC_clust) <- levels(cluster_annots$AUC_clust)

#look at specific clusters
#MSC-like population
chosen <- "MSC_like"
interesting.wrs.msc <- markers.sce.cluster[[chosen]]
interesting.wrs.msc[1:10, 1:4]

chosen <- "Ciliated"
interesting.wrs.cil <- markers.sce.cluster[[chosen]]
interesting.wrs.cil[1:10, 1:4]

chosen <- "Secretory"
interesting.wrs.sec <- markers.sce.cluster[[chosen]]
interesting.wrs.sec[1:10, 1:4]

chosen <- "Early_Secretory"
interesting.wrs.esec <- markers.sce.cluster[[chosen]]
interesting.wrs.esec[1:10, 1:4]

chosen <- "Branch"
interesting.wrs.branch <- markers.sce.cluster[[chosen]]
interesting.wrs.branch[1:10, 1:4]

chosen <- "Immune"
interesting.wrs.imm <- markers.sce.cluster[[chosen]]
interesting.wrs.imm[1:10, 1:4]

## pick off the marker genes
clust_types <- unique(salmon_quants_hg38$cluster.collapse)
marker_genes_by_cluster <- do.call(rbind, lapply(clust_types, function(x) {
  markers <- markers.sce.cluster[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.001,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:7]))
  auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

marker_genes_by_cluster.nocollapse <- do.call(rbind, lapply(unique(salmon_quants_hg38$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.001,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:10]))
  auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

top_marker_genes_by_cluster <- do.call(rbind, lapply(clust_types, function(x) {
  markers <- markers.sce.cluster[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.001,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:7]))
  auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- markers[order(markers$AUC.ave, decreasing = T),]
  markers <- markers[1:10,]
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

top_marker_genes_by_cluster.nocollapse <- do.call(rbind, lapply(unique(salmon_quants_hg38$cluster), function(x) {
  markers <- markers.sce.cluster.nocollapse[[x]]
  ## keep sig
  markers <- markers[markers$FDR < 0.05,]
  ## subset to upper quartile of AUC
  markers$AUC.ave <- rowMeans(as.matrix(markers[,4:10]))
  auc.thresh <- quantile(markers$AUC.ave, probs = seq(0,1,0.05))["95%"]
  markers <- markers[markers$AUC.ave > auc.thresh,]
  markers <- markers[order(markers$AUC.ave, decreasing = T),]
  if (nrow(markers) >= 10) {
    markers <- markers[1:10,]
  }
  markers <- data.frame(marker_genes = rownames(markers),
                        cluster = x)
  return(markers)
}))

## rename clusters
clust_names_nova_nocollapse <- salmon_quants_hg38$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "Ciliated_arm",
    . == 2 ~ "Secretory_arm",
    . == 3 ~ "MSC_like",
    . == 4 ~ "Ciliated",
    . == 5 ~ "Branch",
    . == 6 ~ "Immune",
    . == 7 ~ "Secretory_late",
    . == 8 ~ "Secretory_early"
  ))

salmon_quants_hg38$cluster.names <- clust_names_nova_nocollapse$cluster_cat

# dedupe
marker_genes_by_cluster <- subset(marker_genes_by_cluster, !duplicated(marker_genes_by_cluster$marker_genes))

## make a heatmap
salmon_quants_hg38.markers.new <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% marker_genes_by_cluster$marker_genes,]
salmon_quants_hg38.markers.top <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% top_marker_genes_by_cluster$marker_genes,]
salmon_quants_hg38.markers.top2 <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% top_marker_genes_by_cluster.nocollapse$marker_genes,]

Heatmap(assay(salmon_quants_hg38.markers, "izar"),
        name = "log2(TPM+1)",
        column_split = colData(salmon_quants_hg38.markers.new)$cluster.collapse,
        show_column_names = FALSE,
        show_row_names = FALSE, col = circlize::colorRamp2(breaks = c(seq(0, 6, length.out = 3)), colors = c("#520C90", "goldenrod", "yellow")))

rownames(salmon_quants_hg38.markers.top) <- rowData(salmon_quants_hg38.markers.top)$symbol
Heatmap(assay(salmon_quants_hg38.markers.top, "izar"),
        name = "log2(TPM+1)",
        column_split = colData(salmon_quants_hg38.markers.new)$cluster.collapse,
        #row_split = colData()
        show_column_names = FALSE,
        show_row_names = TRUE, col = circlize::colorRamp2(breaks = c(seq(0, 6, length.out = 3)), colors = c("#520C90", "goldenrod", "yellow")))

rownames(salmon_quants_hg38.markers.top2) <- rowData(salmon_quants_hg38.markers.top2)$symbol
Heatmap(assay(salmon_quants_hg38.markers.top2, "izar"),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38$cluster.names,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE, col = circlize::colorRamp2(breaks = c(seq(0, 10, length.out = 4)), colors = c("#520C90", "#B807EC", "goldenrod", "yellow")))

Heatmap(assay(salmon_quants_hg38.markers.top2, "izar"),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38$cluster.names,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE, col = getJetColors(circlize.cols = TRUE))


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


## let's look at bulk MSC derived from omentum
## load them in
quant.dirs <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/msc_bulk_ens102",
                         pattern = "_quant", full.names = TRUE)
quants <- list.files(quant.dirs, pattern = "quant.sf", full.names = TRUE)
names(quants) <- list.files("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/msc_bulk_ens102",
                            pattern = "_quant")
salmon_quants_hg38_msc_bulk <- import_plate_txis(quants = quants,
                                        t2g = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.tx2gene.tsv",
                                        gtf = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/Homo_sapiens.GRCh38.102.expanded.gtf",
                                        qc = FALSE)
metadata(salmon_quants_hg38_msc_bulk)$origin <- "bulk"
saveRDS(salmon_quants_hg38_msc_bulk, file = "~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/msc_bulk_ens102/bulk_msc_omentum_grch38_ens102_salmon_raw_se.rds")

## combine objects for plotting
salmon_quants_hg38_msc_bulk$cluster.names <- "MSC"
#bulk_genes <- get_ensembl_genes("102")
all(rownames(salmon_quants_hg38) == rownames(salmon_quants_hg38_msc_bulk))
## TRUE
rowData(salmon_quants_hg38_msc_bulk)$symbol <- rowData(salmon_quants_hg38)$symbol
## strip down the single-cell object
salmon_quants_hg38.sub <- salmon_quants_hg38
colData(salmon_quants_hg38.sub) <- colData(salmon_quants_hg38.sub)[,c("sample",
                                                              "NumGenesExpressed",
                                                              "sizeFactor",
                                                              "cluster.names")]
assays(salmon_quants_hg38.sub) <- list(counts = assay(salmon_quants_hg38.sub, "counts"),
                                       logcounts = assay(salmon_quants_hg38.sub, "logcounts"),
                                       spliced = assay(salmon_quants_hg38.sub, "spliced"),
                                       unspliced = assay(salmon_quants_hg38.sub, "unspliced"),
                                       spliced_tpm = assay(salmon_quants_hg38.sub, "spliced_tpm"),
                                       unspliced_tpm = assay(salmon_quants_hg38.sub, "unspliced_tpm"),
                                       tpm = assay(salmon_quants_hg38.sub, "tpm"))
reducedDims(salmon_quants_hg38.sub) <- NULL
combined_objs <- cbind(salmon_quants_hg38.sub, salmon_quants_hg38_msc_bulk)

## plot marker genes
combined_objs.plot <- combined_objs[rowData(combined_objs)$symbol %in% rownames(salmon_quants_hg38.markers.top2),]
rownames(combined_objs.plot) <- rowData(combined_objs.plot)$symbol
combined_objs.plot <- izar_transform(combined_objs.plot)
Heatmap(assay(combined_objs.plot, "izar"),
        name = "log2(TPM+1)",
        column_split = combined_objs.plot$cluster.names,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE, col = getJetColors(circlize.cols = TRUE))

## do a dirty trick with UMAP

## drag in the izar2020 SS2 results
izar_2020 <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/analysis/Izar_et_al_2020_nature_medicine_raw_smartseq2_log_tpm_sce.rds")

rownames(salmon_quants_hg38.markers) <- rowData(salmon_quants_hg38.markers)$symbol
salmon_quants_hg38.markers <- subset(salmon_quants_hg38.markers, !duplicated(rownames(salmon_quants_hg38.markers)))

izar_2020.fte.markers <- izar_2020[rownames(izar_2020) %in% rownames(salmon_quants_hg38.markers),]

## now subset to those found in the izar data
salmon_quants_hg38.markers <- salmon_quants_hg38.markers[rownames(salmon_quants_hg38.markers) %in% rownames(izar_2020.fte.markers),]

izar_match <- match(rownames(salmon_quants_hg38.markers), rownames(izar_2020.fte.markers))
izar_2020.fte.markers <- izar_2020.fte.markers[izar_match,]

## toss out cluster 9 since it has 2 cells in it
izar_2020.fte.markers <- izar_2020.fte.markers[,!izar_2020.fte.markers$clst %in% "9"]

colData(salmon_quants_hg38.markers)$Tissue <- "FTE"
colData(izar_2020.fte.markers)$Tissue <- "HGSOC ascites"

salmon_quants_hg38.markers.all <- as.data.frame(t(cbind(assay(salmon_quants_hg38.markers, "izar"), assay(izar_2020.fte.markers))))
salmon_quants_hg38.markers.all$Tissue <- c(rep("FTE", 308), rep("HGSOC ascites", 1297))
salmon_quants_hg38.markers.all$Cluster <- c(salmon_quants_hg38.markers$cluster.collapse,
                                            paste0("HGSOC_clust_", izar_2020.fte.markers$clst))
salmon_quants_hg38.markers.all.melt <- reshape2::melt(salmon_quants_hg38.markers.all, id.vars = c("Tissue", "Cluster"))

h1 <- Heatmap(scale(as.matrix(t(salmon_quants_hg38.markers.all))),
        name = "log2(TPM+1)",
        column_split = c(salmon_quants_hg38.markers$cluster.collapse,
                         izar_2020.fte.markers$clst),
        show_column_names = FALSE,
        show_row_names = FALSE, col = circlize::colorRamp2(breaks = c(seq(0, 6, length.out = 3)), colors = c("#520C90", "goldenrod", "yellow")),
        column_title_rot = 90)

draw(h1)

## okay let's look at metabolic embeddings now
hallmark_glycolysis <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/gene_sets/hallmark_glycolysis_human.txt",
                                  header = TRUE)
hallmark_oxphos <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/gene_sets/hallmark_oxphos_human.txt",
                              header = TRUE)
hallmark_fa <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/gene_sets/hallmark_lipid_FA_metab_human.txt",
                          header = TRUE)

## tack on gene symbols to gene names
gene_names <- get_ensembl_genes(version = "102", species = "Homo.sapiens")
gene_names.match <- match(rownames(salmon_quants_hg38),
                          gene_names$gene_id)
gene_names <- gene_names[gene_names.match,]
rowData(salmon_quants_hg38)$symbol <- gene_names$symbol

## pull out tpm for each of these pathways and generate a heatmap
glycolysis_tpm <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% hallmark_glycolysis$HALLMARK_GLYCOLYSIS,]
oxphos_tpm <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% hallmark_oxphos$HALLMARK_OXIDATIVE_PHOSPHORYLATION,]
fa_tpm <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% hallmark_fa$HALLMARK_FATTY_ACID_METABOLISM,]

rowData(glycolysis_tpm)$met.pathway <- "Glycolysis"
rowData(oxphos_tpm)$met.pathway <- "OxPhos"
rowData(fa_tpm)$met.pathway <- "FA"

salmon_quants_hg38.met <- rbind(glycolysis_tpm,
                                oxphos_tpm,
                                fa_tpm)
salmon_quants_hg38.met <- salmon_quants_hg38.met[,salmon_quants_hg38.met$cluster %in% c(3, 4, 7, 8)]

## highly variable genes?
assay(salmon_quants_hg38.met, "scaled_tpm") <- scale(assay(salmon_quants_hg38.met, "izar"),
                                                     center = TRUE, scale = FALSE)

assay(salmon_quants_hg38.met, "ratio_us_tpm") <- assay(salmon_quants_hg38.met, "unspliced_tpm") / (assay(salmon_quants_hg38.met, "unspliced_tpm") + assay(salmon_quants_hg38.met, "spliced_tpm"))
assay(salmon_quants_hg38.met, "ratio_us_tpm")[is.nan(assay(salmon_quants_hg38.met, "ratio_us_tpm"))] <- 0

salmon_quants_hg38.met.sub <- salmon_quants_hg38.met[rowMeans(assay(salmon_quants_hg38.met, "ratio_us_tpm")) > 0,]

library(ComplexHeatmap)

Heatmap(assay(salmon_quants_hg38.met.sub, "ratio_us_tpm"),
        row_split = rowData(salmon_quants_hg38.met.sub)$met.pathway,
        column_split = colData(salmon_quants_hg38.met.sub)$cluster,
        show_column_names = FALSE,
        show_row_names = FALSE, col = circlize::colorRamp2(breaks = c(seq(0, 1, length.out = 2)), colors = c("lightgray", "maroon")))

##Hm... this isn't as clear as I'd like.
##switch over to embedding with metabolism genes


## do velocity with scvelo in python :)
## don't use velocessor per-se
salmon_quants_hg38 <- compute_velocity(salmon_quants_hg38, embed = "TSNE",
                                       scvmode = "dynamical")

## rebuild the damn object...
salmon_quants_hg38.stripped <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings_forscvelo.rds")


### detect doublets
library(BiocSingular)
library(scDblFinder)
dbl.dens <- computeDoubletDensity(salmon_quants_hg38, subset.row=top.hvgs, 
                                  dims=ncol(reducedDim(salmon_quants_hg38, "PCAsub")))
summary(dbl.dens)
head(dbl.dens)
salmon_quants_hg38$doubletScore <- dbl.dens
plotReducedDim(salmon_quants_hg38, colour_by = "doubletScore", dimred = "densMAP")
## it shows that a score of around 2 exists for the MSC-like and branch point
## this isn't unexpected and indistinguishable from real biology - needed ERCC spikes...
## hm, looking at the cytoTRACE paper - I wonder if we have gene expression correlates too
salmon_quants_hg38$log_numgenes <- log10(salmon_quants_hg38$NumGenesExpressed)
plotReducedDim(salmon_quants_hg38, colour_by = "log_numgenes", dimred = "densMAP")
## indeed we have a decreasing number of genes from MSC-like to terminally differentiated
## this is likely real


## combine the data
salmon_quants_hg38.hiseq <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")

## look at common marker genes across patients
salmon_quants_hg38.hiseq.markers.top <- salmon_quants_hg38.hiseq[rowData(salmon_quants_hg38.hiseq)$symbol %in% rownames(salmon_quants_hg38.markers.top),]
rownames(salmon_quants_hg38.hiseq.markers.top) <- rowData(salmon_quants_hg38.hiseq.markers.top)$symbol

## combine
combined_izar_tpm_fte <- cbind(assay(salmon_quants_hg38.markers.top, "izar"),
                               assay(salmon_quants_hg38.hiseq.markers.top, "izar"))

combined_izar_tpm_fte.meta <- data.frame()

Heatmap(combined_izar_tpm_fte,
        name = "log2(TPM+1)",
        show_column_names = FALSE,
        show_row_names = TRUE,
        col = circlize::colorRamp2(breaks = c(seq(0, 6, length.out = 3)), colors = c("#520C90", "goldenrod", "yellow")))


## what happens if we look at things with singleR
library(SingleR)
library(celldex)
hpca.se <- celldex::HumanPrimaryCellAtlasData()

salmon_quants_hg38 <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")

to_annot <- assay(salmon_quants_hg38, "logcounts")
rownames(to_annot) <- rowData(salmon_quants_hg38)$symbol
salmon_quants_hg38.annotated <- SingleR(test = to_annot,
                                        ref = hpca.se,
                                        labels = hpca.se$label.fine,
                                        clusters = salmon_quants_hg38$cluster)




colData(salmon_quants_hg38)$ass1 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "ASS1",]
colData(salmon_quants_hg38)$ar <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "AR",]
colData(salmon_quants_hg38)$muc <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MUC1",]
colData(salmon_quants_hg38)$fgfr2 <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FGFR2",]
colData(salmon_quants_hg38)$msln <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "MSLN",]
colData(salmon_quants_hg38)$fkbpl <- assays(salmon_quants_hg38)$izar[rowData(salmon_quants_hg38)$symbol == "FKBPL",]

library(cowplot)
to_plot <- as.data.frame(reducedDim(salmon_quants_hg38, "densMAP"))
colnames(to_plot) <- c("Dim_1", "Dim_2")

ass1 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$ass1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("ASS1 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

pax8 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$pax8)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("PAX8 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ar <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$ar)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("AR gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

muc1 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$muc)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("MUC1 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

fgfr2 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$fgfr2)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("FGFR2 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

foxj1 <- ggplot(to_plot, aes(x = Dim_1, y = Dim_2, color = colData(salmon_quants_hg38)$foxj1)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("FOXJ1 gene expression") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(pax8, foxj1, ass1, ar, muc1, fgfr2, nrow=3, ncol = 2)
cell_markers <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/supp_table_1_marker_genes.txt")


plotReducedDim(salmon_quants_hg38, colour_by = "cluster", dimred = "densMAP")
plotReducedDim(salmon_quants_hg38, colour_by = "ovgp1", dimred = "densMAP")

# forget this... going to do a joint distribution of expression
#BiocManager::install("Nebulosa")
library(Nebulosa)
rownames(salmon_quants_hg38) <- uniquifyFeatureNames(rownames(salmon_quants_hg38),
                                                     rowData(salmon_quants_hg38)$symbol)

top_block <- c("EFNB2", "FGFR2", "NRG3", "VTCN1", "MUC1", "EMX2", "AR",
               "ZNF8", "THSD4", "SLC7A11", "PODXL", "LINC00937", "SNHG3",
               "DHCR24", "ANO1", "FBXO21", "PKHD1L1", "TTYH1",
               "CRTAC1", "OVGP1")

bottom_block <- c("KIRREL1", "GJA1", "TM4SF1", "PAX2", "PAX8", "KRT8",
                  "MSLN", "ITPR3", "ASS1", "C3", "BTNL9", "PAK3",
                  "FAM107A", "ELF3", "EPPK1", "FBXO2")

plot_density(salmon_quants_hg38, features = "PAX8", reduction = "densMAP", joint = TRUE, combine = FALSE)

# not super helpful

# make a heatmap of these blocks of genes and look at overlap with cluster
library(ComplexHeatmap)


salmon_quants_hg38.markers.blocks <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% c(top_block, bottom_block),]

# subset to just secretory and ciliated clusters
salmon_quants_hg38.markers.blocks <- salmon_quants_hg38.markers.blocks[,salmon_quants_hg38.markers.blocks$cluster %in% c(8,7)]

Heatmap(scale(assay(salmon_quants_hg38.markers.blocks, "izar")),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38.markers.blocks$cluster,
        #row_split = salmon_quants_hg38$cluster.names
        show_column_names = FALSE,
        show_row_names = TRUE)


## this is all very confusing - look at the marker genes again
salmon_quants_hg38 <- readRDS("~/Documents/manuscripts/metabolic_flow_fte_2020/salmon_requant/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")

## pull in the clusters that correspond to secretory and ciliated
secretory_c5 <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/marker_expression_tables/cluster_5_marker_exp_table.txt")
secretory_c7 <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/marker_expression_tables/cluster_7_marker_exp_table.txt")
ciliated_c4 <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/marker_expression_tables/cluster_4_marker_exp_table.txt")
ciliated_c6 <- read.delim("~/Documents/manuscripts/metabolic_flow_fte_2020/marker_expression_tables/cluster_6_marker_exp_table.txt")

## subset
secretory_c5.sig.auc <- secretory_c5[secretory_c5$FDR < 0.05,]
secretory_c5.sig.auc <- secretory_c5.sig.auc[secretory_c5.sig.auc$summary.AUC > 0.9,]
secretory_c7.sig.auc <- secretory_c7[secretory_c7$FDR < 0.05,]
secretory_c7.sig.auc <- secretory_c7.sig.auc[secretory_c7.sig.auc$summary.AUC > 0.9,]

## find overlap
intersect(rownames(secretory_c7.sig.auc), rownames(secretory_c5.sig.auc))
# [1] "SNORA73B" "MSLN"     "SNHG3"    "CRISP3"   "C3"       "SNORA81"  "SNORA73A" "SCARNA2"  "PAK3"    
# [10] "PKHD1L1"  "MUC1"     "FBXO21"   "EMX2"     "ASS1"     "FAM107A"  "SNORD46"  "THSD4"    "OVGP1"   
# [19] "AR"       "PODXL"    "ITPR3"    "BTNL9"    "RSPO1"    "FMOD"     "DANT2"    "SNORD33"  "KRT7"    
# [28] "SNHG25"   "PLCB1"    "SNORA24"  "FOSB"     "NFKBIZ"   "SLC38A2"  "FGFR2"    "KRT18"    "CLU"     
# [37] "KIRREL1"  "DHCR24"   "ELF3"     "ANO1"     "KRT8"     "TM4SF1"   "TTYH1"    "PAX2"     "CRTAC1"  
# [46] "SLC7A11"  "GJA1"     "EDN1"     "EPPK1"    "NRG3"     "VTCN1"    "EFNB2"    "SNORD55"  "ADAMTS1" 
# [55] "MUC16"    "MAP3K21"  "CDH12"    "ZNF83"    "C6orf132" "TNKS1BP1"

# look at the mean TPM values across clusters in single patient
# 7 and 8
sec_cells.sce <- salmon_quants_hg38[,salmon_quants_hg38$cluster %in% c(7,8,4)]
sec_cells.sce.markers <- sec_cells.sce[rownames(sec_cells.sce) %in% union(rownames(secretory_c7.sig.auc), rownames(secretory_c5.sig.auc)),]

# threshold expression
sec_cells.sce.markers <- sec_cells.sce.markers[rowMeans(assay(sec_cells.sce.markers, "izar")[,sec_cells.sce.markers$cluster %in% c(7,8)]) > 1,]

# delta between expression in ciliated vs secretory
sec_markers <- rownames(sec_cells.sce.markers)[(rowMeans(assay(sec_cells.sce.markers, "izar")[,sec_cells.sce.markers$cluster %in% c(7,8)]) / 
rowMeans(assay(sec_cells.sce.markers, "izar")[,sec_cells.sce.markers$cluster %in% c(4)])) > 5]

sec_markers <- c(sec_markers, "PAX8", "LINC00937", "RSPO1")
write.table(sec_markers, file = "secretory_markers.txt", quote = F, row.names = F, col.names = F)

## subset
ciliated_c4.sig.auc <- ciliated_c4[ciliated_c4$FDR < 0.05,]
ciliated_c4.sig.auc <- ciliated_c4.sig.auc[ciliated_c4.sig.auc$summary.AUC > 0.95,]
ciliated_c6.sig.auc <- ciliated_c6[ciliated_c6$FDR < 0.05,]
ciliated_c6.sig.auc <- ciliated_c6.sig.auc[ciliated_c6.sig.auc$summary.AUC > 0.95,]

## find overlap
intersect(rownames(ciliated_c6.sig.auc), rownames(ciliated_c4.sig.auc))
# [1] "VWA3B"     "PIFO"      "ZBBX"      "NEK5"      "UBXN10"    "DNAH9"     "CDHR3"     "CFAP54"    "LRRC46"   
# [10] "FHAD1"     "LRRIQ1"    "CFAP70"    "RSPH4A"    "HYDIN"     "DNAH2"     "CFAP45"    "ANKFN1"    "ERICH3"   
# [19] "DNAH6"     "TEKT2"     "TTC29"     "TPPP3"     "DNAH10"    "DLEC1"     "CFAP91"    "CFAP43"    "FAM216B"  
# [28] "WDR49"     "RP1"       "SNTN"      "DRC3"      "MAP3K19"   "MT-ND5"    "C20orf85"  "DNAI1"     "TMEM190"  
# [37] "CCDC17"    "EFCAB12"   "CFAP65"    "TEKT1"     "DTHD1"     "RSPH1"     "CFAP157"   "DNAH5"     "DNAH3"    
# [46] "DNAH12"    "TMEM67"    "CFAP44"    "CFAP52"    "TTC25"     "CCDC39"    "EFCAB1"    "CEP126"    "MNS1"     
# [55] "CASC2"     "DRC1"      "CAPS"      "TMEM231"   "C1orf194"  "CETN2"     "DNAI3"     "DNAI4"     "CCDC40"   
# [64] "CFAP73"    "DNAH11"    "ADGB"      "DNAAF1"    "CFAP46"    "ODF3B"     "CCDC81"    "C6orf118"  "CDH1"     
# [73] "DSP"       "LRRFIP1"   "CCBE1"     "HSPA12A"   "TUB"       "SPTBN1"    "LINC01320"         

## same selection criteria
# look at the mean TPM values across clusters in single patient
# 7 and 8
cil_cells.sce <- salmon_quants_hg38[,salmon_quants_hg38$cluster %in% c(1,4,7,8)]
cil_cells.sce.markers <- cil_cells.sce[rownames(cil_cells.sce) %in% intersect(rownames(ciliated_c4.sig.auc), rownames(ciliated_c6.sig.auc)),]

# threshold expression
cil_cells.sce.markers <- cil_cells.sce.markers[rowMeans(assay(cil_cells.sce.markers, "izar")[,cil_cells.sce.markers$cluster %in% c(1,4)]) > 1,]

# delta between expression in ciliated vs cilretory
cil_markers <- rownames(cil_cells.sce.markers)[(rowMeans(assay(cil_cells.sce.markers, "izar")[,cil_cells.sce.markers$cluster %in% c(1,4)]) / 
                                                  rowMeans(assay(cil_cells.sce.markers, "izar")[,cil_cells.sce.markers$cluster %in% c(7,8)])) > 7]

cil_markers <- c(cil_markers, "FOXJ1")
write.table(cil_markers, file = "ciliated_markers.txt", quote = F, row.names = F, col.names = F)

plot_density(salmon_quants_hg38, features = "PAK3", reduction = "densMAP", joint = TRUE, combine = FALSE)

salmon_quants_hg38.markers.blocks <- salmon_quants_hg38[rowData(salmon_quants_hg38)$symbol %in% c(sec_markers, cil_markers),]

# subset to just secretory and ciliated clusters
#salmon_quants_hg38.markers.blocks <- salmon_quants_hg38.markers.blocks[,salmon_quants_hg38.markers.blocks$cluster %in% c(8,7)]

Heatmap(scale(assay(salmon_quants_hg38.markers.blocks, "izar")),
        name = "log2(TPM+1)",
        column_split = salmon_quants_hg38.markers.blocks$cluster,
        show_column_names = FALSE,
        show_row_names = TRUE, getJetColors(circlize.cols = TRUE)) 



## look at ESR1, BRCA1, PCNA, and MKI67
library(Nebulosa)
rownames(salmon_quants_hg38) <- uniquifyFeatureNames(rownames(salmon_quants_hg38),
                                                     rowData(salmon_quants_hg38)$symbol)

## exclude immune cells
salmon_quants_hg38.noimmune <- salmon_quants_hg38[,!salmon_quants_hg38$cluster %in% 6]

my_theme <- theme(
  plot.title = element_text(hjust = 0.5)
)
esr1 <- plot_density(salmon_quants_hg38.noimmune, features = "ESR1",
                     reduction = "densMAP", combine = FALSE,
                     slot = "izar")
esr1 + my_theme + geom_point(size = 2)

brca1 <- plot_density(salmon_quants_hg38.noimmune, features = "BRCA1",
                     reduction = "densMAP", combine = FALSE,
                     slot = "izar")
brca1 + my_theme + geom_point(size = 2)

foxj1 <- plot_density(salmon_quants_hg38.noimmune, features = "FOXJ1",
                      reduction = "densMAP", combine = FALSE,
                      slot = "izar")
foxj1 + my_theme + geom_point(size = 2)

pax8 <- plot_density(salmon_quants_hg38.noimmune, features = "PAX8",
                      reduction = "densMAP", combine = FALSE,
                      slot = "izar", method = "wkde")
pax8 + my_theme + geom_point(size = 2)

cd45 <- plot_density(salmon_quants_hg38, features = "PTPRC",
                     reduction = "densMAP", combine = FALSE,
                     slot = "izar")
cd45 + my_theme + geom_point(size = 2)

fkbpl <- plot_density(salmon_quants_hg38.noimmune, features = "FKBPL",
                     reduction = "densMAP", combine = FALSE,
                     slot = "logcounts")
fkbpl + my_theme + geom_point(size = 2)

pcna <- plot_density(salmon_quants_hg38.noimmune, features = "PCNA",
                     reduction = "densMAP", combine = FALSE)
pcna + my_theme + geom_point(size = 2)

mki67 <- plot_density(salmon_quants_hg38.noimmune, features = "MKI67",
                     reduction = "densMAP", combine = FALSE)
mki67 + my_theme + geom_point(size = 2)

