# ingest and process the FTE TE results

# patient 1
# mapped using Ensembl and repeatmasker with Ayush's TE pipeline
te_cpm_filt <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/count_matrix_filtered_cpm_TE.tsv.gz",
                          header = TRUE)
rownames(te_cpm_filt) <- te_cpm_filt$Geneid

# patient 2
# mapped using Ensembl and repeatmasker with Ayush's TE pipeline
te_cpm_filt_pat2 <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/pat2/pat2_hiseq_fte_te_cpm_sce_ens101_raw.rds")


# quick filter all zero entries
table(rowSums(te_cpm_filt[,c(1:306)]) > 0)
# FALSE    TRUE 
# 1789937 2748141 

# when it's done with stranded quants...
# FALSE    TRUE
# 2760153 1777925
# nearly inverted...

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
saveRDS(te_cpm_filt.sce, file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/pat1_novaseq_fte_te_cpm_sce_ens101_raw_stranded.rds")
# te_cpm_filt.sce <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/pat1_novaseq_fte_te_cpm_sce_ens101_raw.rds")
te_cpm_filt.sce <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/pat1_novaseq_fte_te_cpm_sce_ens101_ten_perc_filt_line_sine_ltr.rds")
# te_cpm_filt.sce <- readRDS("pat1_novaseq_fte_te_cpm_sce_ens101_ten_perc_filt_line_sine_ltr.intergenic_only.rds")

# now we need to drag in the filtered patient 1 expression object
# so that we can add clusters and filtered cells
pat1_sce <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_novasseq_grch38_salmon_norm_se_with_embeddings.rds")
colnames(pat1_sce) <- gsub("_quant", "", colnames(pat1_sce))

pat2_sce <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")
colnames(pat2_sce) <- gsub("_R1_001.merged_quant|_R1_001_quant", "", colnames(pat2_sce))

# filter to cells found in the TE sce
pat1_sce.filt_vec <- colnames(pat1_sce) %in% colnames(te_cpm_filt.sce)
pat1_sce <- pat1_sce[,pat1_sce.filt_vec]

pat2_sce.filt_vec <- colnames(te_cpm_filt_pat2) %in% colnames(pat2_sce)
te_cpm_filt_pat2 <- te_cpm_filt_pat2[,pat2_sce.filt_vec]

# do it again for the flip
# this whole exercise hurts on a laptop.... do this on HPC!
# start above with the 10% filt sce
# te_cpm_filt.sce.filt_vec <- colnames(te_cpm_filt.sce) %in% colnames(pat1_sce)
# te_cpm_filt.sce <- te_cpm_filt.sce[,te_cpm_filt.sce.filt_vec]

# match and filter cells
pat1_sce.m <- match(colnames(te_cpm_filt.sce),
                    colnames(pat1_sce))
pat1_sce <- pat1_sce[,pat1_sce.m]

# careful
all(colnames(te_cpm_filt.sce) == colnames(pat1_sce))
# TRUE

# tack on the cluster IDs from the pat1_sce
# te_cpm_filt.sce$cluster <- pat1_sce$cluster

# subset to just LINEs, LTRs, and SINES
# te_cpm_filt.sce <- te_cpm_filt.sce[rowData(te_cpm_filt.sce)$repClass %in% c("LINE",
#                                                                             "LTR",
#                                                                             "SINE"),]

# filter TEs to those that are expressed in at least 10% of cells
# keep_genes <- rowMeans(logcounts(te_cpm_filt.sce) > 0) >= 0.1
# FALSE    TRUE
# 2402144   67483

# after doing stranded quants
# FALSE    TRUE
# 1581735   30633

# filter
# te_cpm_filt.sce <- te_cpm_filt.sce[keep_genes, ]

# let's do some NMF on genes and TEs
# first need to filter the genes to those that are expressed in the smallest
# cluster of cells
table(pat1_sce$cluster)
# 1  2  3  4  5  6  7  8 
# 66 53 33 37 32 14 50 20

# let's just call it 30 cells or 10% of the total cells
keep.genes <- rowSums(assay(pat1_sce, "izar") > 0) >= 30
pat1_sce.gene_filt <- pat1_sce[keep.genes,]

# pull off the hvgs
library(scran)

dec <- scran::modelGeneVar(pat1_sce.gene_filt)
top.hvgs <- scran::getTopHVGs(dec, n=3000)

pat1_sce.gene_filt.hvg <- pat1_sce.gene_filt[rownames(pat1_sce.gene_filt) %in% top.hvgs,]

dec.te <- scran::modelGeneVar(te_cpm_filt.sce)
# plot(dec.te$mean, dec.te$total, xlab="Mean log-expression", ylab="Variance")
# curve(metadata(dec.te)$trend(x), col="blue", add=TRUE)
top.te.hvgs <- scran::getTopHVGs(dec.te, n=5000)

te_cpm_filt.sce.hvg <- te_cpm_filt.sce[rownames(te_cpm_filt.sce) %in% top.te.hvgs,]

# run WGCNA on the gene modules
library(WGCNA)
enableWGCNAThreads(nThreads = 40)

genes <- t(assay(pat1_sce.gene_filt.hvg, "izar"))
tes <- t(assay(te_cpm_filt.sce.hvg, "logcounts"))

# seed
set.seed(42)
# choose soft‐thresholding powers
powers <- c(1:10, seq(12,20,2))

sft_genes <- pickSoftThreshold(genes,
                               powerVector = powers,
                               networkType = "signed",
                               corFnc = "bicor",
                               corOptions = list(use = "pairwise.complete.obs"))

sft_tes <- pickSoftThreshold(tes,
                             powerVector = powers,
                             networkType = "signed",
                             corFnc = "bicor",
                             corOptions = list(use = "pairwise.complete.obs"))

power_genes <- sft_genes$powerEstimate # 6
power_tes <- sft_tes$powerEstimate # 8

# try blockwise
bw_genes <- blockwiseModules(
  datExpr        = genes,
  power          = power_genes,
  networkType    = "signed",
  minModuleSize  = 30,
  deepSplit      = 4,
  mergeCutHeight = 0.1,
  pamRespectsDendro = FALSE,
  verbose        = 3
)

colors_genes <- labels2colors(bw_genes$colors)

# colors_genes
# blue     brown     green       red turquoise    yellow
# 362        85      1024       240       407       882

# adjacency & TOM, then modules
# tom_genes <- TOMsimilarityFromExpr(genes,
#                                    power = power_genes,
#                                    networkType = "signed",
#                                    corType = "bicor",
#                                    nThreads = 40)
# dissTOM_g <- 1 - tom_genes
# 
# geneTree <- hclust(as.dist(dissTOM_g), method="average")
# modules_genes <- cutreeDynamic(
#   dendro = geneTree,
#   distM = dissTOM_g,
#   deepSplit = 2,
#   pamRespectsDendro = FALSE,
#   minClusterSize = 30
# )

# ..cutHeight not given, setting it to 0.986  ===>  99% of the (truncated) height range in dendro.

# colors_genes <- labels2colors(modules_genes)
# colors_genes
# black      blue     brown     green       red turquoise    yellow
# 89       722       437       361       174      2838       379

# split out turquoise
# get the indices of turquoise genes
# turq.idx <- which(colors_genes == "turquoise")
# re‐cluster that subset
# tom_gene_sub   <- TOMsimilarityFromExpr(genes[,turq.idx],
#                                         power = power_genes)
# diss_gene_sub  <- 1 - tom_gene_sub
# subGeneTree   <- hclust(as.dist(diss_gene_sub), method="average")
# subGeneMods   <- cutreeDynamic(
#   dendro          = subGeneTree,
#   distM           = diss_gene_sub,
#   deepSplit       = 4,
#   minClusterSize  = 30,
#   cutHeight = 0.9
# )
# colors_sub_genes <- labels2colors(subGeneMods)
# quick qc plot
plotDendroAndColors(
  dendro     = geneTree,
  colors     = colors_genes,
  groupLabels= "Dynamic Modules",
  main       = "Dendrogram + module colors"
)
# 6 modules
# looks okay, let's move on

# try blockwise
bw_tes <- blockwiseModules(
  datExpr          = tes,
  power             = power_tes,
  networkType       = "signed",
  minModuleSize     = 30,
  deepSplit         = 2,
  mergeCutHeight    = 0.1,
  pamRespectsDendro = FALSE,
  verbose           = 3
)

colors_tes <- labels2colors(bw_tes$colors)

# colors_tes
# black      blue     brown     green   magenta      pink       red turquoise
# 55       190       177      3385       165       845        45        50
# yellow
# 88

save(bw_genes, bw_tes,
     pat1_sce.gene_filt.hvg,
     te_cpm_filt.sce.hvg,
     file = "wgcna_tes_genes_pat1_blockwise_use_this.rda")
# load("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/wgcna_tes_genes_pat1_blockwise_use_this.rda")

# tom_tes  <- TOMsimilarityFromExpr(tes,
#                                   power = power_tes,
#                                   networkType = "signed",
#                                   corType = "bicor",
#                                   nThreads = 40)
# dissTOM_t <- 1 - tom_tes
# 
# teTree  <- hclust(as.dist(dissTOM_t), method="average")
# modules_tes <- cutreeDynamic(
#   dendro      = teTree,
#   distM       = dissTOM_t,
#   deepSplit   = 2,
#   pamRespectsDendro=FALSE,
#   minClusterSize   = 30
# )

# ..cutHeight not given, setting it to 0.995  ===>  99% of the (truncated) height range in dendro.

# colors_tes <- labels2colors(modules_tes)

# colors_tes
# blue     brown     green       red turquoise    yellow
# 313       308       207        42      1915       215

# split out turquoise
# get the indices of turquoise genes
# turq.idx <- which(colors_tes == "turquoise")
# re‐cluster that subset
# tom_te_sub   <- TOMsimilarityFromExpr(tes[,turq.idx],
#                                         power = power_tes)
# diss_te_sub  <- 1 - tom_te_sub
# subTETree   <- hclust(as.dist(diss_te_sub), method="average")
# subTEMods   <- cutreeDynamic(
#   dendro          = subTETree,
#   distM           = diss_te_sub,
#   deepSplit       = 4,
#   minClusterSize  = 30
# )
# colors_sub_tes <- labels2colors(subTEMods)
# 4) module eigengenes
# ME_genes <- moduleEigengenes(genes,
#                              colors = colors_genes,
#                              softPower = power_genes)$eigengenes
# ME_tes   <- moduleEigengenes(tes,
#                              colors = colors_tes,
#                              softPower = power_tes)$eigengenes

ME_genes <- bw_genes$MEs
ME_tes <- bw_tes$MEs

# project to hiseq data
common_genes <- intersect(rownames(pat2_sce),
                          names(bw_genes$colors))
common_tes <- intersect(rownames(te_cpm_filt_pat2),
                        names(bw_tes$colors))

# subset
pat2_sce.sub <- pat2_sce[rownames(pat2_sce) %in% common_genes,]
te_cpm_filt_pat2.sub <- te_cpm_filt_pat2[rownames(te_cpm_filt_pat2) %in% common_tes,]

cols_te <- bw_tes$colors[names(bw_tes$colors) %in% common_tes]

## pull off genes and tes
genes_hiseq <- t(assay(pat2_sce.sub, "izar"))
tes_hiseq <- t(assay(te_cpm_filt_pat2.sub, "logcounts"))

ME_gene_hiseq <- moduleEigengenes(genes_hiseq, colors = bw_genes$colors, 
                                  align = "along average", 
                                  excludeGrey = TRUE)
ME_te_hiseq <- moduleEigengenes(tes_hiseq, colors = cols_te, 
                                align = "along average", 
                                excludeGrey = TRUE)
save.image("~/Documents/manuscripts/storm_seq/te_quant/fte/pat2/projected_me_genes_image_20250916.rda")

# Optional: sign-align to a template ME (e.g., reference cell-type means)
# for each module m:
# if cor(ME_new[,m], ME_template[,m], use="p") < 0 -> ME_new[,m] <- -ME_new[,m]
# 5) cross‐correlate eigengenes
#    You may want to match sample‐order (rows) across both sets!
stopifnot(rownames(ME_genes) == rownames(ME_tes))
stopifnot(rownames(ME_gene_hiseq) == rownames(ME_te_hiseq))

cor_mat <- bicor(ME_genes, ME_tes,
                 use = "pairwise.complete.obs")
p_mat   <- corPvalueStudent(cor_mat,
                            nSamples = nrow(ME_genes))

cor_mat2 <- bicor(ME_gene_hiseq$eigengenes, ME_te_hiseq$eigengenes,
                 use = "pairwise.complete.obs")
p_mat2   <- corPvalueStudent(cor_mat2,
                            nSamples = nrow(ME_gene_hiseq$eigengenes))

# flatten
p_vec <- as.vector(p_mat)
names(p_vec) <- paste(
  rep(colnames(ME_genes), each = ncol(ME_tes)),
  rep(colnames(ME_tes), times = ncol(ME_genes)),
  sep = "_vs_"
)

p_vec2 <- as.vector(p_mat2)
names(p_vec2) <- paste(
  rep(colnames(ME_gene_hiseq$eigengenes), each = ncol(ME_te_hiseq$eigengenes)),
  rep(colnames(ME_te_hiseq$eigengenes), times = ncol(ME_gene_hiseq$eigengenes)),
  sep = "_vs_"
)

# adjust
fdr_vec <- p.adjust(p_vec, method = "BH")
fdr_vec2 <- p.adjust(p_vec2, method = "BH")

# reconstruct a matrix
fdr_mat <- matrix(fdr_vec,
                  nrow = ncol(ME_genes),
                  ncol = ncol(ME_tes),
                  dimnames = list(colnames(ME_genes),
                                  colnames(ME_tes)))

fdr_mat2 <- matrix(fdr_vec2,
                  nrow = ncol(ME_gene_hiseq$eigengenes),
                  ncol = ncol(ME_te_hiseq$eigengenes),
                  dimnames = list(colnames(ME_gene_hiseq$eigengenes),
                                  colnames(ME_te_hiseq$eigengenes)))

# name the genes and TEs again
rowData(te_cpm_filt.sce.hvg)$wgcna_colors <- bw_tes$colors
rowData(pat1_sce.gene_filt.hvg)$wgcna_colors <- bw_genes$colors

rowData(te_cpm_filt_pat2.sub)$wgcna_colors <- ME_te_hiseq$validColors
rowData(pat2_sce.sub)$wgcna_colors <- ME_gene_hiseq$validColors

# cluster the eigengenes
# metree_tes <- hclust(as.dist(1 - bicor(ME_tes)), method="average")
# plot(metree_tes, main="TE module ME clustering")
# mergeCutHeight <- 0.5
# abline(h = mergeCutHeight, col="red", lty=2)

# merge
# mergeRes_tes <- mergeCloseModules(
#   tes,
#   colors_tes,
#   cutHeight = mergeCutHeight,
#   verbose   = 3
# )
# colors_tes_merged <- mergeRes_tes$colors
# length(unique(colors_tes_merged))
# table(colors_tes_merged)

# 6) visualize as a heatmap
# filter out the grey modules
cor_mat <- cor_mat[!rownames(cor_mat) %in% "MEgrey",
                   !colnames(cor_mat) %in% "MEgrey"]
fdr_mat <- fdr_mat[!rownames(fdr_mat) %in% "MEgrey",
                   !colnames(fdr_mat) %in% "MEgrey"]
ME_tes <- ME_tes[,!colnames(ME_tes) %in% "MEgrey"]
ME_genes <- ME_genes[,!colnames(ME_genes) %in% "MEgrey"]


labeledHeatmap(
  Matrix = cor_mat,
  xLabels = colnames(ME_tes),
  yLabels = colnames(ME_genes),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = paste(signif(cor_mat,2), "\n(", signif(fdr_mat,1), ")", sep=""),
  zlim = c(-1,1),
  main = "Correlations: Gene Modules vs TE Modules"
)

labeledHeatmap(
  Matrix = cor_mat2,
  xLabels = colnames(ME_te_hiseq$eigengenes),
  yLabels = colnames(ME_gene_hiseq$eigengenes),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = paste(signif(cor_mat2,2), "\n(", signif(fdr_mat2,1), ")", sep=""),
  zlim = c(-1,1),
  main = "Correlations: Gene Modules vs TE Modules"
)

# rescale the module eigengenes for joint density plotting
ME_genes.rescale <- as.data.frame(lapply(ME_genes, function(x) {
  scales::rescale(x, to = c(0,1))
}),
row.names = rownames(ME_genes))

ME_tes.rescale <- as.data.frame(lapply(ME_tes, function(x) {
  scales::rescale(x, to = c(0,1))
}),
row.names = rownames(ME_tes))

# pat 2
ME_genes.rescale2 <- as.data.frame(lapply(ME_gene_hiseq$eigengenes, function(x) {
  scales::rescale(x, to = c(0,1))
}),
row.names = rownames(ME_gene_hiseq$eigengenes))

ME_tes.rescale2 <- as.data.frame(lapply(ME_te_hiseq$eigengenes, function(x) {
  scales::rescale(x, to = c(0,1))
}),
row.names = rownames(ME_te_hiseq$eigengenes))

library(scran)
library(scater)
reducedDim(te_cpm_filt.sce.hvg, "densMAP") <- reducedDim(pat1_sce.gene_filt.hvg,
                                                         "densMAP")
colData(te_cpm_filt.sce.hvg) <- cbind(colData(te_cpm_filt.sce.hvg),
                                      ME_tes.rescale)
colData(pat1_sce.gene_filt.hvg) <- cbind(colData(pat1_sce.gene_filt.hvg),
                                         ME_genes.rescale)

reducedDim(te_cpm_filt_pat2.sub, "PCA") <- reducedDim(pat2_sce.sub,
                                                          "PCA")
colData(te_cpm_filt_pat2.sub) <- cbind(colData(te_cpm_filt_pat2.sub),
                                      ME_tes.rescale2)
colData(pat2_sce.sub) <- cbind(colData(pat2_sce.sub),
                                         ME_genes.rescale2)

# YES!
# we have module eigengenes composed of TEs that correspond
# to different cell types
# testing...
colData(pat1_sce)$notch <- assays(pat1_sce)$logcounts[rowData(pat1_sce)$symbol == "NOTCH1",]
colData(pat2_sce)$notch <- assays(pat2_sce)$logcounts[rowData(pat2_sce)$symbol == "NOTCH1",]
colData(pat1_sce)$olfm4 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "OLFM4",]
colData(pat1_sce)$znf382 <- assays(pat1_sce)$logcounts[rowData(pat1_sce)$symbol == "ZNF382",]
colData(pat1_sce)$foxa2 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "FOXA2",]
colData(pat1_sce)$pkd1l1 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "PKD1L1",]
colData(pat1_sce)$foxj1 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "FOXJ1",]

notch <- plotReducedDim(pat1_sce, dimred = "densMAP", colour_by = "notch")
znf <- plotReducedDim(pat1_sce, dimred = "densMAP", colour_by = "znf382")
foxj1 <- plotReducedDim(pat1_sce, dimred = "densMAP", colour_by = "foxj1")
foxa2 <- plotReducedDim(pat1_sce, dimred = "densMAP", colour_by = "foxa2")

# mitos
colData(pat1_sce)$nrf1 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "NRF1",]
colData(pat1_sce)$nrf2 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "NFE2L2",]
colData(pat1_sce)$tfam <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "TFAM",]
colData(pat1_sce)$pgc1a <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "PPARGC1A",]

colData(pat1_sce)$mtrnr1 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "MT-RNR1",]
colData(pat1_sce)$mtnd1 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "MT-ND1",]
colData(pat1_sce)$mtnd4 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "MT-ND4",]
colData(pat1_sce)$mtnd5 <- assays(pat1_sce)$izar[rowData(pat1_sce)$symbol == "MT-ND5",]

# mito joint
mito_genes <- rowData(pat1_sce)$symbol[grep("^MT-", rowData(pat1_sce)$symbol)]
colData(pat1_sce)$mito_exp <- colSums(assays(pat1_sce)$izar[rowData(pat1_sce)$symbol %in% mito_genes,])

# Split: mt mRNA (OXPHOS) vs rRNA/tRNA
mt_mrna <- c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6",
             "MT-CO3","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-CYB","MT-ND6")
mt_rrna_trna <- setdiff(mito_genes, mt_mrna)

sym <- rowData(pat1_sce)$symbol
sym2 <- rowData(pat2_sce)$symbol
ix_all  <- which(sym %in% mito_genes)
ix_all2  <- which(sym2 %in% mito_genes)
ix_mrna <- which(sym %in% mt_mrna)
ix_mrna2 <- which(sym2 %in% mt_mrna)
ix_non  <- setdiff(seq_along(sym), ix_all)
ix_non2  <- setdiff(seq_along(sym2), ix_all2)

# If you already have TPM:
tpm <- assays(pat1_sce)$counts
tpm2 <- assays(pat2_sce)$counts

# If you have counts and gene lengths (kb) in rowData(sce)$length_kb:
# tpm <- t(t(t(assay(sce, "counts")) / rowData(sce)$length_kb) * 1e6 / 
#          colSums(t(assay(sce, "counts") / rowData(sce)$length_kb)))

# Per-cell totals
tot_tpm    <- colSums(tpm)
tot_tpm2    <- colSums(tpm2)
mt_all_tpm <- colSums(tpm[ix_all, , drop=FALSE])
mt_all_tpm2 <- colSums(tpm2[ix_all2, , drop=FALSE])
mt_mrna_tpm<- colSums(tpm[ix_mrna, , drop=FALSE])
mt_mrna_tpm2<- colSums(tpm2[ix_mrna2, , drop=FALSE])
mt_rrna_tpm<- colSums(tpm[match(mt_rrna_trna, sym), , drop=FALSE], na.rm=TRUE)
mt_rrna_tpm2<- colSums(tpm2[match(mt_rrna_trna, sym2), , drop=FALSE], na.rm=TRUE)

frac_mt_all  <- mt_all_tpm / tot_tpm
frac_mt_all2  <- mt_all_tpm2 / tot_tpm2
frac_mt_mrna <- mt_mrna_tpm / tot_tpm
frac_mt_mrna2 <- mt_mrna_tpm2 / tot_tpm2
frac_mt_rrna <- mt_rrna_tpm / tot_tpm
frac_mt_rrna2 <- mt_rrna_tpm2 / tot_tpm2

# Log2 TPM for selected genes (for density plots)
genes_of_interest <- mt_mrna  # or mt_all
log2_tpm <- tpm[match(genes_of_interest, sym), , drop=FALSE]
log2_tpm2 <- tpm2[match(genes_of_interest, sym2), , drop=FALSE]

# Store in colData
colData(pat1_sce)$frac_mt_all  <- frac_mt_all
colData(pat1_sce)$frac_mt_mrna <- frac_mt_mrna
colData(pat1_sce)$frac_mt_rrna <- frac_mt_rrna

colData(pat2_sce)$frac_mt_all  <- frac_mt_all2
colData(pat2_sce)$frac_mt_mrna <- frac_mt_mrna2
colData(pat2_sce)$frac_mt_rrna <- frac_mt_rrna2

# make a plot stratified by cell type and color by mt fraction
pat1_sce$cluster.collapse <- pat1_sce$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "Ciliated",
    . == 2 ~ "Secretory",
    . == 3 ~ "UCFP",
    . == 4 ~ "Ciliated",
    . == 5 ~ "Branch",
    . == 6 ~ "Immune",
    . == 7 ~ "Secretory",
    . == 8 ~ "Secretory"
  ))
pat1_sce$cluster.collapse <- pat1_sce$cluster.collapse$cluster_cat

pat2_sce$cluster.collapse <- pat2_sce$cluster %>%
  as_tibble() %>%
  mutate(cluster_cat = case_when(
    . == 1 ~ "UCFP",
    . == 2 ~ "Secretory",
    . == 3 ~ "Ciliated",
    . == 4 ~ "Secretory",
    . == 5 ~ "Secretory",
    . == 6 ~ "Ciliated"
  ))
pat2_sce$cluster.collapse <- pat2_sce$cluster.collapse$cluster_cat

# let's call cell cycle too - see if the higher mito content goes with cell cycle
# No, this fails because of the modeling approach
# we should do the dye cycle violet experiment during sorting for scRNA prep


to_plot_mt <- data.frame(mt_frac = c(pat1_sce$frac_mt_all,
                                     pat2_sce$frac_mt_all),
                         cluster = c(pat1_sce$cluster.collapse,
                                     pat2_sce$cluster.collapse))

ggplot(to_plot_mt, aes(x = cluster, y = mt_frac,
                       color = mt_frac, group = factor(cluster))) +
  geom_jitter(width = 0.25, alpha = 0.6) +
  stat_summary(geom = "crossbar", fun = "median", linewidth = 1) +
  scale_color_viridis_c(option = "plasma",
                        name = "Read Frac") +
  geom_hline(yintercept = c(0.05, 0.10), linetype = "dashed", linewidth = 1.2) +
  #geom_hline(yintercept = 0.10, linetype = "dashed", linewidth = 1) +
  ylim(c(0,0.25)) +
  ylab("Fraction mitochondrial reads") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

plotReducedDim(pat1_sce, dimred = "densMAP", colour_by = "frac_mt_mrna")

library(Nebulosa)
# priming for mitochondrial biogenesis
pgc1a <- plot_density(pat1_sce, "pgc1a") +
  ggtitle("PGC1a Expression Density") +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

# consolidation with ciliagenesis tfs foxj1 and foxa2
nrf_joint <- plot_density(pat1_sce, c("nrf1","nrf2"), joint = TRUE)[[3]] +
  ggtitle("NRF1/2 Expression Density") +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

foxj1 <- plot_density(pat1_sce, "foxj1") +
  ggtitle("FOXJ1 Expression Density") +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

foxa2 <- plot_density(pat1_sce, "foxa2") +
  ggtitle("FOXA2 Expression Density") +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

plot_grid(pgc1a, nrf_joint, foxj1, foxa2, nrow = 1, ncol = 4, align = "hv")

plotReducedDim(te_cpm_filt.sce.hvg, dimred = "densMAP", colour_by = "MEturquoise")
plotReducedDim(pat1_sce.gene_filt.hvg, dimred = "densMAP", colour_by = "notch")

plotReducedDim(te_cpm_filt_pat2.sub, dimred = "PCA", colour_by = "MEturquoise")
plotReducedDim(pat2_sce.sub, dimred = "PCA", colour_by = "MEturquoise")

# do the plots to show joint density

# need to combine objects
ME_tes.annot <- paste0(colnames(ME_tes.rescale), "_TE")
ME_tes.annot.mat <- ME_tes.rescale
colnames(ME_tes.annot.mat) <- ME_tes.annot
colData(pat1_sce.gene_filt.hvg) <- cbind(colData(pat1_sce.gene_filt.hvg),
                                         ME_tes.annot.mat)

ME_tes.annot2 <- paste0(colnames(ME_tes.rescale2), "_TE")
ME_tes.annot.mat2 <- ME_tes.rescale2
colnames(ME_tes.annot.mat2) <- ME_tes.annot2
colData(pat2_sce.sub) <- cbind(colData(pat2_sce.sub),
                                         ME_tes.annot.mat2)

to_plot <- data.frame(dim_1 = reducedDim(pat1_sce.gene_filt.hvg, "densMAP")[,1],
                      dim_2 = reducedDim(pat1_sce.gene_filt.hvg, "densMAP")[,2],
                      ciliated = pat1_sce.gene_filt.hvg$MEturquoise *
                        pat1_sce.gene_filt.hvg$MEturquoise_TE,
                      secretory = pat1_sce.gene_filt.hvg$MEyellow *
                        pat1_sce.gene_filt.hvg$MEbrown_TE,
                      secretory_sub = pat1_sce.gene_filt.hvg$MEgreen *
                        pat1_sce.gene_filt.hvg$MEyellow_TE,
                      ucfp = pat1_sce.gene_filt.hvg$MEbrown *
                        pat1_sce.gene_filt.hvg$MEblack_TE,
                      immune = pat1_sce.gene_filt.hvg$MEblue *
                        pat1_sce.gene_filt.hvg$MEred_TE)

to_plot <- data.frame(dim_1 = reducedDim(pat2_sce.sub, "PCA")[,1],
                      dim_2 = reducedDim(pat2_sce.sub, "PCA")[,2],
                      ciliated = pat2_sce.sub$MEturquoise *
                        pat2_sce.sub$MEturquoise_TE,
                      secretory = pat2_sce.sub$MEyellow *
                        pat2_sce.sub$MEbrown_TE,
                      secretory_sub = pat2_sce.sub$MEgreen *
                        pat2_sce.sub$MEyellow_TE,
                      ucfp = pat2_sce.sub$MEbrown *
                        pat2_sce.sub$MEblack_TE,
                      immune = pat2_sce.sub$MEblue *
                        pat2_sce.sub$MEred_TE)

library(ggplot2)
library(cowplot)

# ciliated: ME_te - turquoise; ME_gene - turquoise
ciliated_plot <- ggplot(to_plot, aes(x = dim_1, y = dim_2, color = ciliated)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Joint Ciliated Gene and TE Eigengenes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# secretory: ME_te - brown; ME_gene - yellow
secretory_plot <- ggplot(to_plot, aes(x = dim_1, y = dim_2, color = secretory)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Joint Secretory Gene and TE Eigengenes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# top secretory: ME_te - yellow; ME_gene - green
sub_secretory_plot <- ggplot(to_plot, aes(x = dim_1, y = dim_2, color = secretory_sub)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Joint 'Upper' Secretory Gene and TE Eigengenes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# UCFP
ucfp_plot <- ggplot(to_plot, aes(x = dim_1, y = dim_2, color = ucfp)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Joint UCFP Gene and TE Eigengenes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# Immune
immune_plot <- ggplot(to_plot, aes(x = dim_1, y = dim_2, color = immune)) +
  geom_point() +
  scale_color_viridis_c(begin = 0.2, name = "Log2 Expression") +
  theme_half_open() +
  ylab("Dimension 2") +
  xlab("Dimension 1") +
  #coord_cartesian(xlim = c(-8, 8), ylim = c(-16, 16)) +
  ggtitle("Joint Immune Gene and TE Eigengenes") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

plot_grid(ciliated_plot, secretory_plot, ucfp_plot,
          align = "hv", nrow = 1, ncol = 3)

# now we need to perform within cell type bicorrelation to try and see if these actually
# co-vary or are largely "markers" for a given cell type
## rename clusters
# clust_names_nova_nocollapse <- salmon_quants_hg38$cluster %>%
#   as_tibble() %>%
#   mutate(cluster_cat = case_when(
#     . == 1 ~ "Ciliated_arm",
#     . == 2 ~ "Secretory_arm",
#     . == 3 ~ "MSC_like",
#     . == 4 ~ "Ciliated",
#     . == 5 ~ "Branch",
#     . == 6 ~ "Immune",
#     . == 7 ~ "Secretory_late",
#     . == 8 ~ "Secretory_early"
#   ))
# ciliated clusters 1 and 4
ciliated_cells <- pat1_sce.gene_filt.hvg[,pat1_sce.gene_filt.hvg$cluster %in% c(1, 4)]
bicorAndPvalue(ciliated_cells$MEturquoise,
      ciliated_cells$MEturquoise_TE,
      use = "pairwise.complete.obs")
# 0.5103966

# ciliated clusters 3 and 6 pat 2
ciliated_cells2 <- pat2_sce.sub[,pat2_sce.sub$cluster %in% c(3, 6)]
bicorAndPvalue(ciliated_cells2$MEturquoise,
               ciliated_cells2$MEturquoise_TE,
               use = "pairwise.complete.obs")
# 0.5103966

# secretory clusters 2, 7, and 8
secretory_cells <- pat1_sce.gene_filt.hvg[,pat1_sce.gene_filt.hvg$cluster %in% c(2, 7, 8)]
bicorAndPvalue(secretory_cells$MEyellow,
      secretory_cells$MEbrown_TE,
      use = "pairwise.complete.obs")
# 0.2504924

# what about just the upper secretory?
secretory_cells <- pat1_sce.gene_filt.hvg[,pat1_sce.gene_filt.hvg$cluster %in% c(2, 7, 8)]
bicor(secretory_cells$MEgreen,
      secretory_cells$MEyellow_TE,
      use = "pairwise.complete.obs")
# 0.6564936

# what about p53 to just the secretory cells?
sec.p53 <- secretory_cells[rownames(secretory_cells) %in% "ENSG00000141510",]
bicor(as.vector(assay(sec.p53, "izar")),
      secretory_cells$MEgreen,
      use = "pairwise.complete.obs")
# 0.09332631 nope!

bicor(as.vector(assay(sec.p53, "izar")),
      secretory_cells$MEyellow,
      use = "pairwise.complete.obs")
# 0.1907646 not a hub gene in the yellow module WITHIN secretory cells...

# what about just the ucfp?
ucfp_cells <- pat1_sce.gene_filt.hvg[,pat1_sce.gene_filt.hvg$cluster %in% c(3)]
bicorAndPvalue(ucfp_cells$MEbrown,
      ucfp_cells$MEblack_TE,
      use = "pairwise.complete.obs")
# 0.2563214

# prioritize the ciliated and upper secretory modules to look for dynamic TE regulation
# with the corresponding gene modules

# correlate these with pseudotime
load("~/Documents/manuscripts/storm_seq/fte_analysis/principle_curves_pat1_pat2_with_pseudotime_20250616.rda")

# correct the names
rownames(to_plot.nova) <- gsub("_quant",
                               "",
                               rownames(to_plot.nova))
to_plot.nova <- to_plot.nova[rownames(to_plot.nova) %in%
                               rownames(ME_tes),]
ME_tes.sub <- ME_tes.rescale[rownames(ME_tes.rescale) %in% rownames(to_plot.nova),]
to_plot.nova.m <- match(rownames(ME_tes.sub),
                        rownames(to_plot.nova))
to_plot.nova <- to_plot.nova[to_plot.nova.m,]
all(rownames(to_plot.nova) == rownames(ME_tes.sub))
# TRUE
ME_genes.sub <- ME_genes.rescale[rownames(ME_genes.rescale) %in% rownames(to_plot.nova),]


# correlate with pseudotime
te_pseudotime <- bicorAndPvalue(ME_tes.sub,
                                to_plot.nova$pseudotime,
                                use = "pairwise.complete.obs")
te_pseudotime.df <- data.frame(
  module = rownames(te_pseudotime$bicor),
  bicor   = as.numeric(te_pseudotime$bicor),
  p_val = as.numeric(te_pseudotime$p),
  q_val = as.numeric(p.adjust(as.numeric(te_pseudotime$p),
                              method = "BH"))
)

gene_pseudotime <- bicorAndPvalue(ME_genes.sub,
                                  to_plot.nova$pseudotime,
                                  use = "pairwise.complete.obs")
gene_pseudotime.df <- data.frame(
  module = rownames(gene_pseudotime$bicor),
  bicor   = as.numeric(gene_pseudotime$bicor),
  p_val = as.numeric(gene_pseudotime$p),
  q_val = as.numeric(p.adjust(as.numeric(gene_pseudotime$p),
                              method = "BH"))
)

# note that immune cells were _not_ included in the pseudotime
# inference, so it's expected that a gene or TE program would be
# anti-correlated with pseudotime (MEblue [gene] and MEred [TE])

# what about these trends
# upper secretory is _strong_

# order the pseudotime
to_plot.nova.order <- to_plot.nova[order(to_plot.nova$pseudotime),]
ME_genes.sub.order <- ME_genes.sub[order(to_plot.nova$pseudotime),] 
ME_tes.sub.order <- ME_tes.sub[order(to_plot.nova$pseudotime),] 

# now only keep the secretory from ucfp 
to_plot.nova.order.sec <- to_plot.nova.order[to_plot.nova.order$clusters %in% c(3, 5, 2, 7, 8),]
ME_genes.sub.order.sec <- ME_genes.sub.order[to_plot.nova.order$clusters %in% c(3, 5, 2, 7, 8),] 
ME_tes.sub.order.sec <- ME_tes.sub.order[to_plot.nova.order$clusters %in% c(3, 5, 2, 7, 8),] 

ME_genes.sub.order.sec$MEgreen <- scales::rescale(ME_genes.sub.order.sec$MEgreen)
ME_tes.sub.order.sec$MEyellow <- scales::rescale(ME_tes.sub.order.sec$MEyellow)
# write.table(to_plot.nova.order.sec, file = "~/Downloads/pseudotime_secretory.txt", quote = F, row.names = F, col.names = T, sep = '\t')
# write.table(ME_genes.sub.order.sec, file = "~/Downloads/gene_modules_secretory.txt", quote = F, row.names = F, col.names = T, sep = '\t')
# write.table(ME_tes.sub.order.sec, file = "~/Downloads/te_modules_secretory.txt", quote = F, row.names = F, col.names = T, sep = '\t')
# now only keep the ciliated from ucfp 
# to_plot.nova.order.cil <- to_plot.nova.order[to_plot.nova.order$clusters %in% c(3, 5, 1, 4),]
# ME_genes.sub.order.cil <- ME_genes.sub.order[to_plot.nova.order$clusters %in% c(3, 5, 1, 4),] 
# ME_tes.sub.order.cil <- ME_tes.sub.order[to_plot.nova.order$clusters %in% c(3, 5, 1, 4),] 

to_plot.nova.order.cil <- to_plot.nova.order[to_plot.nova.order$clusters %in% c(5, 1, 4),]
ME_genes.sub.order.cil <- ME_genes.sub.order[to_plot.nova.order$clusters %in% c(5, 1, 4),] 
ME_tes.sub.order.cil <- ME_tes.sub.order[to_plot.nova.order$clusters %in% c(5, 1, 4),] 

pat1_sce.gene_filt.hvg.turq <- logcounts(pat1_sce.gene_filt.hvg)[rownames(pat1_sce.gene_filt.hvg) %in% names(bw_genes$colors)[bw_genes$colors %in% "turquoise"],]
te_cpm_filt.sce.hvg.turq <- logcounts(te_cpm_filt.sce.hvg)[rownames(te_cpm_filt.sce.hvg) %in% names(bw_tes$colors)[bw_tes$colors %in% "turquoise"],]

# get the cells we are after
pat1_sce.gene_filt.hvg.turq <- pat1_sce.gene_filt.hvg.turq[,colnames(pat1_sce.gene_filt.hvg.turq) %in% rownames(ME_genes.sub.order.cil)]
te_cpm_filt.sce.hvg.turq <- te_cpm_filt.sce.hvg.turq[,colnames(te_cpm_filt.sce.hvg.turq) %in% rownames(ME_tes.sub.order.cil)]

# reorder to pseudotime
gene.m <- match(rownames(to_plot.nova.order.cil),
                colnames(pat1_sce.gene_filt.hvg.turq))
pat1_sce.gene_filt.hvg.turq <- pat1_sce.gene_filt.hvg.turq[,gene.m]

te.m <- match(rownames(to_plot.nova.order.cil),
              colnames(te_cpm_filt.sce.hvg.turq))
te_cpm_filt.sce.hvg.turq <- te_cpm_filt.sce.hvg.turq[,te.m]

# side step to look at module TEs that change over pseudotime
## ---------- (A) Module-level trajectories ----------
# Per-cell module means (average across features)
te_cell_mean   <- colMeans(te_cpm_filt.sce.hvg.turq)
gene_cell_mean <- colMeans(pat1_sce.gene_filt.hvg.turq)

library(tidyr)

df_cell <- tibble(
  pseudotime = pt,
  TE         = te_cell_mean,
  Gene       = gene_cell_mean
) |>
  pivot_longer(c(TE, Gene), names_to = "type", values_to = "expr")

# Option A2: quantile-binned means (guarantees non-empty bins)
nbins <- 15
qbreaks <- quantile(pt, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
cuts <- cut(pt, breaks = unique(qbreaks), include.lowest = TRUE)   # remove dup breakpoints if pt not diverse

bin_center <- tapply(pt, cuts, mean)
te_bin     <- tapply(te_cell_mean, cuts, mean)
gene_bin   <- tapply(gene_cell_mean, cuts, mean)

df_bins <- tibble(
  pseudo_x = as.numeric(bin_center),
  TE       = as.numeric(te_bin),
  Gene     = as.numeric(gene_bin)
) |>
  pivot_longer(c(TE, Gene), names_to = "type", values_to = "mean_expr") |>
  drop_na()

p_bins <- ggplot(df_bins, aes(pseudo_x, mean_expr, color = type)) +
  geom_line(size = 1.1) +
  geom_point(size = 1.4) +
  labs(title = "Module mean expression (quantile-binned)",
       x = "Pseudotime (bin center)", y = "Mean expression (module features)", color = "Feature") +
  theme_minimal(base_size = 12)
p_bins

## ---------- (B) Per-feature trend classification ----------
# Fit a simple trend for each feature vs pseudotime and BH-correct p-values.
# You can swap LM for Spearman if you prefer rank-robust trends.

fit_trend_mat <- function(M, label, pt, spearman_exact = FALSE) {
  # M: features x cells (columns already ordered by pseudotime)
  # pt: numeric pseudotime (length = ncol(M))
  purrr::map_dfr(seq_len(nrow(M)), function(i){
    y <- as.numeric(M[i, ])
    
    ## Linear trend (for reference)
    m  <- lm(y ~ pt)
    co <- summary(m)$coefficients
    lm_slope <- unname(co["pt","Estimate"])
    lm_pval  <- unname(co["pt","Pr(>|t|)"])
    
    ## Spearman correlation (rank-based, robust to monotone nonlinearities)
    # exact=FALSE avoids very slow exact p-values with many ties
    ct <- suppressWarnings(cor.test(pt, y, method = "spearman", exact = spearman_exact))
    rho <- unname(ct$estimate)
    sp_pval <- unname(ct$p.value)
    
    tibble::tibble(
      feature   = rownames(M)[i] %||% paste0(label, "_", i),
      type      = label,
      lm_slope  = lm_slope,
      lm_pval   = lm_pval,
      sp_rho    = rho,
      sp_pval   = sp_pval
    )
  })
}

library(purrr)
## --- Run for TEs and genes ---
trend_te   <- fit_trend_mat(te_cpm_filt.sce.hvg.turq,   "TE",   pt)
trend_gene <- fit_trend_mat(pat1_sce.gene_filt.hvg.turq, "Gene", pt)
trend_tbl  <- dplyr::bind_rows(trend_te, trend_gene) |>
  dplyr::mutate(
    lm_padj = p.adjust(lm_pval, method = "BH"),
    sp_padj = p.adjust(sp_pval, method = "BH")
  )

## --- Choose classification metric ---
# Option A: classify by Spearman (recommended for pseudotime trends)
trend_by_spearman <- trend_tbl |>
  dplyr::mutate(class = dplyr::case_when(
    sp_padj < 0.05 & sp_rho >  0 ~ "activated",
    sp_padj < 0.05 & sp_rho <  0 ~ "repressed",
    TRUE                         ~ "flat"
  ))

# Option B: classify by LM slope (kept for reference)
trend_by_lm <- trend_tbl |>
  dplyr::mutate(class = dplyr::case_when(
    lm_padj < 0.05 & lm_slope >  0 ~ "activated",
    lm_padj < 0.05 & lm_slope <  0 ~ "repressed",
    TRUE                           ~ "flat"
  ))

## --- Summary counts (Spearman-based) ---
trend_by_spearman |>
  dplyr::count(type, class, name = "n") |>
  dplyr::arrange(type, dplyr::desc(n)) |>
  print()

trend_by_lm |>
  dplyr::count(type, class, name = "n") |>
  dplyr::arrange(type, dplyr::desc(n)) |>
  print()

## --- Volcano-style plots ---
# Spearman version
ggplot(trend_by_spearman, aes(sp_rho, -log10(sp_padj), color = class, shape = type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.4) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.3) +
  geom_point(alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "Feature trends vs pseudotime (Spearman)",
       x = "Spearman rho", y = "-log10(FDR)", color = "Class", shape = "Type")

# LM version (optional)
ggplot(trend_by_lm, aes(lm_slope, -log10(lm_padj), color = class, shape = type)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.4) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.3) +
  geom_point(alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "Feature trends vs pseudotime (LM)",
       x = "Slope (LM)", y = "-log10(FDR)", color = "Class", shape = "Type")

library(ppcor)   # partial correlation

# 1) keep only significant activated/repressed TEs
sig_tes <- trend_by_spearman %>%
  dplyr::filter(type == "TE", sp_padj < 0.05, class %in% c("activated","repressed")) %>%
  dplyr::select(feature, class)

# 2) subset matrix to those TEs (and check names)
keep_features <- intersect(sig_tes$feature, rownames(te_cpm_filt.sce.hvg.turq))
expr_sub <- te_cpm_filt.sce.hvg.turq[keep_features, , drop = FALSE]

# 3) matrix -> long tibble: one row per (cell, TE)
df_te <- as_tibble(t(expr_sub), rownames = "cell") %>%   # now cells x TEs
  tidyr::pivot_longer(
    cols = -cell,
    names_to = "feature",
    values_to = "expr"
  ) %>%
  # add pseudotime (assumes colnames(expr_te_mat) == cell names)
  dplyr::mutate(pseudotime = pt[match(cell, colnames(te_cpm_filt.sce.hvg.turq))]) %>%
  # add class (activated/repressed)
  dplyr::left_join(sig_tes, by = "feature") %>%
  # keep ordered by pseudotime for plotting
  dplyr::arrange(pseudotime)

# Quantile breaks (guarantees non-empty bins unless ties collapse levels)
nbins <- 15
qbr   <- quantile(df_te$pseudotime, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
qbr   <- unique(qbr)  # remove duplicate breakpoints if pseudotime has ties

df_te_b <- df_te %>%
  mutate(bin = cut(pseudotime, breaks = qbr, include.lowest = TRUE)) %>%
  filter(!is.na(bin))

# 1) Bin centers computed ONCE
bin_centers <- df_te_b %>%
  group_by(bin) %>%
  summarise(bin_pt = mean(pseudotime), .groups = "drop")

# 2) Mean per TE within each bin
te_bin_means <- df_te_b %>%
  group_by(class, feature, bin) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

# 3) Class-level mean ± sd per bin, then add bin centers
df_bins <- te_bin_means %>%
  group_by(class, bin) %>%
  summarise(
    bin_mean = mean(mean_expr),
    bin_sd   = sd(mean_expr),
    .groups  = "drop"
  ) %>%
  left_join(bin_centers, by = "bin") %>%
  arrange(class, bin)

# Plot
ggplot(df_bins, aes(bin_pt, bin_mean, color = class, fill = class)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = bin_mean - bin_sd, ymax = bin_mean + bin_sd),
              alpha = 0.20, color = NA) +
  theme_minimal(base_size = 12) +
  labs(title = "Mean TE trajectories (quantile-binned, warning-free)",
       x = "Pseudotime (bin center)", y = "Mean TE expression")

# Spaghetti with per-TE fitted lines + class average overlay
df_te_repress <- df_te[df_te$class %in% "repressed",]
ggplot() +
  # raw per-cell expression, colored by class only
  geom_point(data = df_te,
             aes(x = pseudotime, y = expr, color = class),
             alpha = 0.15, size = 0.5) +
  # per-TE fitted lines (light, grouped by feature)
  geom_smooth(data = df_te,
              aes(x = pseudotime, y = expr, group = feature, color = class),
              method = "gam", formula = y ~ s(x, k = 3),
              se = FALSE, alpha = 0.2, linewidth = 0.6) +
  # optional class-level overlay (bold mean trajectory)
  # geom_smooth(data = df_te,
  #             aes(x = pseudotime, y = expr, color = class),
  #             method = "gam", formula = y ~ s(x, k = 3),
  #             se = FALSE, linewidth = 1.4) +
  facet_wrap(~ class, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "TE trajectories with raw per-cell expression + per-TE fits",
       x = "Pseudotime", y = "TE expression (logcounts)")

# monotonic decrease along trajectory of genes and a _sharp_ increase for the TEs near early secretory transition
plot(to_plot.nova.order.sec$pseudotime, ME_genes.sub.order.sec[,"MEgreen"],   type="l", col="black", ylim=c(0,1))
lines(to_plot.nova.order.sec$pseudotime, ME_tes.sub.order.sec[  ,"MEyellow"],col="firebrick")
legend("topright", legend=c("Genes","TEs"), col=c("black","firebrick"), lwd=2)

plot(to_plot.nova.order.cil$pseudotime, ME_genes.sub.order.cil[,"MEturquoise"],   type="l", col="black", ylim=c(0,1))
lines(to_plot.nova.order.cil$pseudotime, ME_tes.sub.order.cil[  ,"MEturquoise"],col="firebrick")
legend("topright", legend=c("Genes","TEs"), col=c("black","firebrick"), lwd=2)


# quickly check if p53 is bicorrelated with this TE module
genes.p53 <- genes[,colnames(genes) %in% "ENSG00000141510"]
genes.dusp5 <- genes[,colnames(genes) %in% "ENSG00000138166"]
genes.pax8 <- genes[,colnames(genes) %in% "ENSG00000125618"]

genes.foxj1 <- genes[,colnames(genes) %in% "ENSG00000129654"]
genes.ccdc17 <- genes[,colnames(genes) %in% "ENSG00000159588"]
genes.ccdc78 <- genes[,colnames(genes) %in% "ENSG00000162004"]

# quick sanity check for bicor with the yellow module in which it's found
bicorAndPvalue(genes.p53,
               ME_genes$MEyellow,
               use = "pairwise.complete.obs")
# bicor: 0.4529265
# pval: 7.821952e-17
# sanity has been checked :)

genes.p53 <- data.frame(p53_exp = genes.p53,
                        dusp5_exp = genes.dusp5,
                        pax8_exp = genes.pax8)

genes.p53 <- genes.p53[rownames(genes.p53) %in% rownames(ME_genes.sub.order.sec),]
genes.p53.m <- match(rownames(ME_genes.sub.order.sec),
                     rownames(genes.p53))
genes.p53 <- genes.p53[genes.p53.m,]

bicorAndPvalue(genes.p53$p53_exp,
               ME_genes.sub.order.sec$MEgreen,
               use = "pairwise.complete.obs")

genes.p53$sample <- rownames(genes.p53)
genes.p53$p53_exp <- scales::rescale(genes.p53$p53_exp, to = c(0,1))
genes.p53$dusp5_exp <- scales::rescale(genes.p53$dusp5_exp, to = c(0,1))
genes.p53$pax8_exp <- scales::rescale(genes.p53$pax8_exp, to = c(0,1))

# ciliated genes
genes.cil <- data.frame(foxj1_exp = genes.foxj1,
                        ccdc17_exp = genes.ccdc17,
                        ccdc78_exp = genes.ccdc78)

genes.cil <- genes.cil[rownames(genes.cil) %in% rownames(ME_genes.sub.order.cil),]
genes.cil.m <- match(rownames(ME_genes.sub.order.cil),
                     rownames(genes.cil))
genes.cil <- genes.cil[genes.cil.m,]

bicorAndPvalue(genes.cil$foxj1_exp,
               ME_genes.sub.order.cil$MEturquoise,
               use = "pairwise.complete.obs")

genes.cil$sample <- rownames(genes.cil)
genes.cil$foxj1_exp <- scales::rescale(genes.cil$foxj1_exp, to = c(0,1))
genes.cil$ccdc17_exp <- scales::rescale(genes.cil$ccdc17_exp, to = c(0,1))
genes.cil$ccdc78_exp <- scales::rescale(genes.cil$ccdc78_exp, to = c(0,1))

# 2) Select what we need and join by sample
library(dplyr)
library(tidyverse)
sec_df <- ME_genes.sub.order.sec %>%
  select(MEgreen) %>%
  rownames_to_column("sample") %>%
  left_join(
    ME_tes.sub.order.sec %>% select(MEyellow) %>% rownames_to_column("sample"),
    by = "sample"
  ) %>%
  left_join(
    to_plot.nova.order.sec %>% select(pseudotime, clusters) %>% rownames_to_column("sample"),
    by = "sample"
  ) %>%
  # left_join(
  #   genes.p53 %>% select(p53_exp, dusp5_exp, pax8_exp) %>% rownames_to_column("sample"),
  #   by = "sample"
  # )
  left_join(
    genes.p53 %>% select(p53_exp, pax8_exp) %>% rownames_to_column("sample"),
    by = "sample"
  )

# 4) Pivot to long format
# sec_df_long <- sec_df %>%
#   select(sample, pseudotime, clusters, MEgreen, MEyellow,
#          p53_exp, dusp5_exp, pax8_exp) %>%
#   pivot_longer(
#     cols      = c(MEgreen, MEyellow, p53_exp,
#                   dusp5_exp, pax8_exp),
#     names_to  = "module",
#     values_to = "score"
#   ) %>%
#   mutate(
#     module = recode(module,
#                     MEgreen  = "Genes (MEgreen)",
#                     MEyellow = "TEs (MEyellow)",
#                     p53_exp = "p53 Expression",
#                     dusp5_exp = "DUSP5 Expression",
#                     pax8_exp = "PAX8 Expression")
#   )

sec_df_long <- sec_df %>%
  select(sample, pseudotime, clusters, MEgreen, MEyellow,
         p53_exp, pax8_exp) %>%
  pivot_longer(
    cols      = c(MEgreen, MEyellow, p53_exp,
                  pax8_exp),
    names_to  = "module",
    values_to = "score"
  ) %>%
  mutate(
    module = recode(module,
                    MEgreen  = "Genes (MEgreen)",
                    MEyellow = "TEs (MEyellow)",
                    p53_exp = "p53 Expression",
                    # dusp5_exp = "DUSP5 Expression",
                    pax8_exp = "PAX8 Expression")
  )

ggplot(sec_df_long, aes(x = pseudotime, y = score, color = module)) +
  geom_line(linewidth = 1, alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3, size = 2) +
  scale_color_manual(values = c("Genes (MEgreen)" = "black",
                                "TEs (MEyellow)" = "red",
                                "p53 Expression" = "darkgray",
                                # "DUSP5 Expression" = "maroon",
                                "PAX8 Expression" = "goldenrod")) +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  theme_bw(12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

cil_df <- ME_genes.sub.order.cil %>%
  select(MEturquoise) %>%
  rownames_to_column("sample") %>%
  left_join(
    ME_tes.sub.order.cil %>% select(MEturquoise) %>% rownames_to_column("sample"),
    by = "sample"
  ) %>%
  left_join(
    to_plot.nova.order.cil %>% select(pseudotime, clusters) %>% rownames_to_column("sample"),
    by = "sample"
  )
  # left_join(
  #   genes.cil %>% select(foxj1_exp) %>% rownames_to_column("sample"),
  #   by = "sample"
  # )

cil_df_long <- cil_df %>%
  select(sample, pseudotime, clusters, MEturquoise.x, MEturquoise.y) %>%
  pivot_longer(
    cols      = c(MEturquoise.x, MEturquoise.y),
                  # foxj1_exp),
                  # ccdc17_exp,
                  # ccdc78_exp),
    names_to  = "module",
    values_to = "score"
  ) %>%
  mutate(
    module = recode(module,
                    MEturquoise.x  = "Genes (MEturquoise)",
                    MEturquoise.y = "TEs (MEturquoise)")
                    #foxj1_exp = "FOXJ1 Expression")
                    # ccdc17_exp = "CCDC17 Expression",
                    # ccdc78_exp = "CCDC78 Expression")
  )

centers    <- c(0.54, 0.67)
half_width <- 0.05

shade_df <- data.frame(
  xmin = pmax(0, centers - half_width),
  xmax = pmin(1, centers + half_width),
  ymin = -Inf, ymax = Inf
)

ggplot(cil_df_long, aes(x = pseudotime, y = score, color = module)) +
  # shaded windows behind everything
  geom_rect(data = shade_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "steelblue", alpha = 0.2) +
  geom_vline(xintercept = centers, linetype = 3, alpha = 0.6) +
  geom_line(linewidth = 1, alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE, span = 0.1, size = 2) +
  scale_color_manual(values = c("Genes (MEturquoise)" = "black",
                                "TEs (MEturquoise)"   = "red")) +
  # coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  theme_bw(12) +
  theme(
    plot.title  = element_text(hjust = 0.5),
    axis.text   = element_text(size = 16),
    axis.title  = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title= element_blank()
  )

# Quantile bins (equal cell counts)
qs <- quantile(cil_df$pseudotime, probs = seq(0, 1, by = 0.1), na.rm = TRUE)

min_n <- 5  # require at least this many cells per bin
loc_df <- do.call(rbind, lapply(seq_len(length(qs) - 1), function(i){
  sel <- cil_df$pseudotime >= qs[i] & cil_df$pseudotime < qs[i+1]
  # include right edge in last bin
  if (i == length(qs) - 1) sel <- cil_df$pseudotime >= qs[i] & cil_df$pseudotime <= qs[i+1]
  n <- sum(sel)
  if (n >= min_n) {
    r <- bicor(cil_df$MEturquoise.y[sel],
               cil_df$MEturquoise.x[sel],
               use = "pairwise.complete.obs")
    p <- corPvalueStudent(r, nSamples = n)
  } else {
    r <- NA_real_; p <- NA_real_
  }
  data.frame(bin=i, start=qs[i], end=qs[i+1],
             mid=mean(cil_df$pseudotime[sel]),
             n=n, r=r, p=p)
}))

# FDR
p.adjust(loc_df$p, method="BH")
# [1] 0.63841486 0.64419569 0.63841486 0.64419569 0.01757208 0.63841486 0.22059126 0.63841486 0.76868640
# [10] 0.03399628

# Sliding-window correlation (smoother than 10 fixed bins)
library(mgcv)
library(WGCNA)

ord  <- order(cil_df$pseudotime)
pt   <- cil_df$pseudotime[ord]
x    <- cil_df$MEturquoise.x[ord]
y    <- cil_df$MEturquoise.y[ord]


# window size and midpoint
k    <- ceiling(length(pt) * 0.1)      # e.g., 10% of cells
half <- floor(k/2)

pt_mid <- sapply(seq_along(pt), function(i){
  i1 <- max(1, i-half); i2 <- min(length(pt), i+half)
  median(pt[i1:i2], na.rm = TRUE)
})

fit_te  <- gam(y ~ s(pt, k=6))   # x, y are standardized
fit_ge  <- gam(x ~ s(pt, k=6))
rx <- residuals(fit_te, type="response")
ry <- residuals(fit_ge, type="response")

# rolling correlation on residuals (centered)
r_roll_resid <- sapply(seq_along(pt), function(i){
  i1 <- max(1, i-half); i2 <- min(length(pt), i+half)
  segx <- rx[i1:i2]; segy <- ry[i1:i2]
  if (sum(is.finite(segx)&is.finite(segy)) < 6) return(NA_real_)
  if (sd(segx,na.rm=TRUE)==0 || sd(segy,na.rm=TRUE)==0) return(NA_real_)
  as.numeric(bicor(segx, segy, use="pairwise.complete.obs"))
})

plot(pt_mid, r_roll_resid, type="l",
     xlab = "Pseudotime (window midpoints)",
     ylab = "Local bicor of residuals (TE ↔ Gene)")
abline(h = 0.5, lty = 2)

# per-index window size actually used (edges are smaller)
n_win <- sapply(seq_along(pt), function(i){
  i1 <- max(1, i-half); i2 <- min(length(pt), i+half); i2 - i1 + 1
})

# p-values and FDR for residual rolling corr
p_roll  <- WGCNA::corPvalueStudent(r_roll_resid, nSamples = n_win)
padj    <- p.adjust(p_roll, method = "BH")

# find local peaks above a threshold (e.g., r >= 0.5) that are FDR<0.05
thr <- 0.5
peaks <- which(diff(sign(diff(r_roll_resid))) == -2) + 1
keep  <- peaks[r_roll_resid[peaks] >= thr & padj[peaks] < 0.05]
data.frame(
  pt = round(pt_mid[keep], 3),
  r  = round(r_roll_resid[keep], 3),
  FDR = signif(padj[keep], 3)
)

# pt     r     FDR
# 1 0.57 0.696 0.03620
# 2 0.69 0.849 0.00271

df <- data.frame(pt=pt_mid, r=r_roll_resid, FDR=padj)
df[df < 0] <- 0.0001
ggplot(df, aes(pt, r)) +
  geom_hline(yintercept=0.5, linetype=2) +
  geom_line() +
  geom_point(data=subset(df, r>=0.5 & FDR<0.05), size = 2, color="firebrick") +
  theme_half_open() +
  ylab("Rolling bicor of residuals (TE ↔ Gene)") +
  xlab("Pseudotime (early to late)") +
  ylim(c(0, 1.0)) +
  #xlim(c(0, 1.0)) +
  ggtitle("Residual bicorrelation of\nciliated cells along pseudotime") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# pseudo median point for plotting
# 0.5759835
# 0.6944302
layer_points <- data.frame(pc1 = c(79.3356, 62.59667),
                           pc2 = c(40.50884, -5.939918))

ggplot(to_plot.nova, aes(x = pc1, y = pc2, color = pseudotime)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "inferno",
                        limits = c(0,1),
                        breaks = seq(0,1,1)) +
  geom_point(data = layer_points, size = 8, color = "black") +
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



safe_window_lag_bicor <- function(xw, yw, L, min_n = 6){
  stopifnot(length(xw) == length(yw))
  n <- length(xw)
  if (n < min_n) return(NA_real_)
  if (L > 0) {                 # genes lead TEs by L
    if (n - L < min_n) return(NA_real_)
    ix <- seq_len(n - L)
    iy <- (L + 1L):n
  } else if (L < 0) {          # TEs lead genes by -L
    L2 <- -L
    if (n - L2 < min_n) return(NA_real_)
    ix <- (L2 + 1L):n
    iy <- seq_len(n - L2)
  } else {                     # L == 0
    ix <- iy <- seq_len(n)
  }
  as.numeric(bicor(xw[ix], yw[iy], use = "pairwise.complete.obs"))
}

# 2) Best local lag within a centered window ----------------------------------
local_best_lag <- function(i, x, y, pt, half, Lfrac = 0.3, min_n = 6){
  i1 <- max(1L, i - half)
  i2 <- min(length(x), i + half)
  xw <- x[i1:i2]
  yw <- y[i1:i2]
  n  <- length(xw)
  if (n < min_n) return(list(rmax = NA_real_, L = NA_integer_, lead_pt = NA_real_, center = median(pt[i1:i2])))
  Lmax <- max(1L, floor(Lfrac * n))
  Ls   <- seq.int(-Lmax, Lmax)
  cors <- sapply(Ls, function(L) safe_window_lag_bicor(xw, yw, L, min_n = min_n))
  best <- which.max(cors)
  step <- median(diff(pt[i1:i2]), na.rm = TRUE)
  list(
    rmax    = cors[best],
    L       = Ls[best],                 # +L: genes lead TEs
    lead_pt = Ls[best] * step,          # lag in pseudotime units
    center  = median(pt[i1:i2], na.rm=TRUE)
  )
}

# 3) Apply at your peak indices (or across grid) -------------------------------
# Inputs you already have:
# pt = ordered pseudotime vector
# rx, ry = residuals from your GAM fits (genes = x, TEs = y)
# half = floor(ceiling(length(pt)*0.10)/2)   # e.g., 10% window
# keep = integer indices of peaks in r_roll_resid

lag_out <- lapply(keep, local_best_lag, x = rx, y = ry, pt = pt, half = half)

lag_df <- data.frame(
  pt_mid = round(vapply(lag_out, `[[`, numeric(1), "center"), 3),
  r0     = round(r_roll_resid[keep], 3),
  rmax   = round(vapply(lag_out, `[[`, numeric(1), "rmax"), 3),
  lag    = vapply(lag_out, `[[`, integer(1), "L"),           # +lag: genes lead
  lead_pt= round(vapply(lag_out, `[[`, numeric(1), "lead_pt"), 3)
)
lag_df

# pull off the genes/TEs that went into each module
table(bw_tes$colors)
# black      blue     brown     green      grey 
# 50       190       177        88      3385 
# pink       red turquoise    yellow 
# 45        55       845       165

table(bw_genes$colors)
# blue     brown     green      grey turquoise 
# 407       362        85       882      1024 
# yellow 
# 240

# split things off for each gene/TE module
# pull in the GTF to annotate genes though
gtf <- rtracklayer::import("~/Documents/manuscripts/storm_seq/rrna/Homo_sapiens.GRCh38.101.ercc92patched.gtf.gz")
int_tes <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/intergenic_intronic_tes.txt.gz",
                      header = TRUE)

# green gene module is interesting as it seems to correlate with
# non-terminally differentiating cells...
green_gene_module <- rownames(pat1_sce.gene_filt.hvg)[rowData(pat1_sce.gene_filt.hvg)$wgcna_colors %in% "green"]
write.table(green_gene_module,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/green_gene_module_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')
# ciliated
turquoise_gene_module <- rownames(pat1_sce.gene_filt.hvg)[rowData(pat1_sce.gene_filt.hvg)$wgcna_colors %in% "turquoise"]
write.table(turquoise_gene_module,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/turquoise_gene_module_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

turuoise_te_module <- names(bw_tes$colors)[bw_tes$colors %in% "turquoise"]
write.table(turuoise_te_module,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/turquoise_te_module_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')
turquoise_te_module.tes <- int_tes[int_tes$name %in% turuoise_te_module$V1,]
write.table(turquoise_te_module.tes[,c("seqnames", "start", "end", "name")],
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/turquoise_te_module_genes.bed",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

# secretory
yellow_gene_module <- rownames(pat1_sce.gene_filt.hvg)[rowData(pat1_sce.gene_filt.hvg)$wgcna_colors %in% "yellow"]
write.table(yellow_gene_module,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/yellow_gene_module_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')


# ucfp
brown_gene_module <- rownames(pat1_sce.gene_filt.hvg)[rowData(pat1_sce.gene_filt.hvg)$wgcna_colors %in% "brown"]
write.table(brown_gene_module,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/brown_gene_module_genes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

## look for TEs in the "bursting" range for TF enrichment
burst_TEs_by_peak <- function(
    te_expr,                     # cells x features OR features x cells
    pseudotime,                  # numeric vector (#cells)
    peak_windows,                # data.frame with columns start, end
    test = c("t","wilcox"),
    alt = c("two.sided","greater","less"),
    min_cells_in = 20,
    min_cells_out = 20,
    assume_log2p1 = TRUE,        # your case
    compute_raw_effects = TRUE,  # back-transform for effect sizes
    eps = 1e-8,
    fdr_method = "BH",
    fdr_by = c("contrast","global"),  # FDR per peak (default) or across all peaks
    min_abs_delta_log = 0,       # optional effect-size gate on log2p1 delta
    min_log2FC_raw   = 0         # optional effect-size gate on raw log2FC
){
  test <- match.arg(test)
  alt  <- match.arg(alt)
  fdr_by <- match.arg(fdr_by)
  
  # Ensure cells x features orientation
  if (nrow(te_expr) != length(pseudotime) && ncol(te_expr) == length(pseudotime)) {
    te_expr <- t(te_expr)
  }
  stopifnot(nrow(te_expr) == length(pseudotime))
  stopifnot(all(c("start","end") %in% colnames(peak_windows)))
  
  # helper
  back_tf <- function(M) if (assume_log2p1) (2^M - 1) else M
  
  # build masks for each peak
  in_list <- lapply(seq_len(nrow(peak_windows)), function(i){
    pseudotime >= peak_windows$start[i] & pseudotime <= peak_windows$end[i]
  })
  
  # choose test
  p_fun <- switch(test,
                  t = function(v, in_idx, out_idx){
                    vin <- v[in_idx]; vout <- v[out_idx]
                    if (sum(is.finite(vin)) < min_cells_in || sum(is.finite(vout)) < min_cells_out) return(NA_real_)
                    suppressWarnings(t.test(vin, vout, alternative = alt)$p.value)
                  },
                  wilcox = function(v, in_idx, out_idx){
                    vin <- v[in_idx]; vout <- v[out_idx]
                    if (sum(is.finite(vin)) < min_cells_in || sum(is.finite(vout)) < min_cells_out) return(NA_real_)
                    suppressWarnings(wilcox.test(vin, vout, alternative = alt)$p.value)
                  }
  )
  
  # run per peak vs REST (which includes other peaks)
  per_peak <- lapply(seq_along(in_list), function(i){
    in_idx  <- in_list[[i]]
    out_idx <- !in_idx & is.finite(pseudotime)
    if (sum(in_idx)  < min_cells_in)  stop(sprintf("Peak %d: too few in-window cells.", i))
    if (sum(out_idx) < min_cells_out) stop(sprintf("Peak %d: too few out-of-window cells.", i))
    
    TE_in  <- te_expr[in_idx, , drop=FALSE]
    TE_out <- te_expr[out_idx,, drop=FALSE]
    
    # log2p1 summaries
    mean_in_log  <- colMeans(TE_in,  na.rm=TRUE)
    mean_out_log <- colMeans(TE_out, na.rm=TRUE)
    delta_log    <- mean_in_log - mean_out_log
    
    # raw summaries (optional)
    if (compute_raw_effects) {
      mean_in_raw  <- colMeans(back_tf(TE_in),  na.rm=TRUE)
      mean_out_raw <- colMeans(back_tf(TE_out), na.rm=TRUE)
      delta_raw    <- mean_in_raw - mean_out_raw
      log2FC_raw   <- log2((mean_in_raw + eps)/(mean_out_raw + eps))
    } else {
      mean_in_raw <- mean_out_raw <- delta_raw <- log2FC_raw <- NA_real_
    }
    
    # prevalence on log2p1 scale
    prop_in  <- colMeans(TE_in  > 0, na.rm=TRUE)
    prop_out <- colMeans(TE_out > 0, na.rm=TRUE)
    
    # p-values
    pvals <- apply(te_expr, 2, p_fun, in_idx=in_idx, out_idx=out_idx)
    
    data.frame(
      feature      = colnames(te_expr),
      peak_id      = i,
      start        = peak_windows$start[i],
      end          = peak_windows$end[i],
      n_in         = sum(in_idx),
      n_out        = sum(out_idx),
      mean_in_log  = mean_in_log,
      mean_out_log = mean_out_log,
      delta_log    = delta_log,
      mean_in_raw  = mean_in_raw,
      mean_out_raw = mean_out_raw,
      delta_raw    = delta_raw,
      log2FC_raw   = log2FC_raw,
      prop_in      = prop_in,
      prop_out     = prop_out,
      pval         = pvals,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  
  res <- do.call(rbind, per_peak)
  
  # FDR
  if (fdr_by == "contrast") {
    res$FDR <- ave(res$pval, res$peak_id, FUN = function(p) p.adjust(p, method = fdr_method))
  } else {
    res$FDR <- p.adjust(res$pval, method = fdr_method)
  }
  
  # flags
  res$hit_two_sided <- res$FDR <= 0.05
  res$hit_up_log    <- res$FDR <= 0.05 & res$delta_log > 0
  res$hit_down_log  <- res$FDR <= 0.05 & res$delta_log < 0
  
  # optional effect-size gating
  if (!is.null(min_abs_delta_log) && min_abs_delta_log > 0) {
    res$hit_two_sided <- res$hit_two_sided & abs(res$delta_log) >= min_abs_delta_log
    res$hit_up_log    <- res$hit_up_log    & abs(res$delta_log) >= min_abs_delta_log
    res$hit_down_log  <- res$hit_down_log  & abs(res$delta_log) >= min_abs_delta_log
  }
  if (!is.null(min_log2FC_raw) && isTRUE(compute_raw_effects) && min_log2FC_raw > 0) {
    res$hit_two_sided <- res$hit_two_sided & abs(res$log2FC_raw) >= min_log2FC_raw
    res$hit_up_log    <- res$hit_up_log    &     res$log2FC_raw  >= min_log2FC_raw
    res$hit_down_log  <- res$hit_down_log  &     res$log2FC_raw  <= -min_log2FC_raw
  }
  
  # order: strongest absolute delta first within each peak
  res <- res[order(res$peak_id, -abs(res$delta_log), res$FDR), ]
  res
}

# pseudotime is a numeric vector (0..1), aligned to rows of the matrices
# pt     r     FDR
# 1 0.57 0.696 0.03620
# 2 0.69 0.849 0.00271
wins <- data.frame(start = c(0.54, 0.66),
                   end   = c(0.60, 0.72))

# TE expression
te_cpm_filt.sce.hvg.cil <- te_cpm_filt.sce.hvg[,te_cpm_filt.sce.hvg$cluster %in% c(1,4)]
te_cpm_filt.sce.hvg.cil <- te_cpm_filt.sce.hvg.cil[rowData(te_cpm_filt.sce.hvg.cil)$wgcna_colors %in% "turquoise",]
te_turq_exp <- t(logcounts(te_cpm_filt.sce.hvg.cil))

detRate <- colMeans(te_turq_exp > 0)               # fraction of cells expressed
baseMean <- colMeans(te_turq_exp)                  # mean level (on the same scale you test)
trendVar <- apply(te_turq_exp, 2, var)             # or MAD across pseudotime-binned means

keep <- detRate >= 0.10 & baseMean >= 0.01 & trendVar > quantile(trendVar, 0.25)
te_expr_f <- te_expr[, keep, drop = FALSE]

te_turq_exp.filt <- te_turq_exp[,detRate >= 0.2913]

# pseudotime
te_cil_pt <- to_plot.nova[rownames(to_plot.nova) %in% rownames(te_turq_exp),]
all(rownames(te_cil_pt) == rownames(te_turq_exp))
# TRUE
# te_expr_turq: cells x TEs (turquoise only); gene_expr_turq not needed for bursting
out <- burst_TEs_by_peak(
  te_expr = te_turq_exp.filt,                  # log2(exp+1)
  pseudotime = te_cil_pt$pseudotime,
  peak_windows = wins,
  test="wilcox", alt="two.sided",
  min_cells_in=10, min_cells_out=10,
  assume_log2p1=TRUE, compute_raw_effects=TRUE,
  fdr_by = "contrast"
)

head(out$main)       # ranked TE "bursters" across the union of peaks
head(out$by_window)  # per-peak-window deltas/log2FCs


# genes         : cells × genes matrix (log-normalized)
# tes           : cells × TEs   matrix (log-normalized)
# ME_genes      : cells × gene‐modules eigengenes (named “ME<color>”)
# ME_tes        : cells × TE‐modules   eigengenes
# colors_genes  : named vector gene → moduleColor
# colors_tes    : named vector TE   → moduleColor

getTopPairs <- function(geneMod, teMod, n = 50) {
  # 0) grab the colors
  colors_genes <- bw_genes$colors
  colors_tes <- bw_tes$colors
  
  # 1) grab the feature names in each module
  genes_interest <- names(colors_genes)[colors_genes %in% geneMod]
  tes_interest   <- names(colors_tes)[colors_tes %in% teMod]
  
  # 1a) filter out constant features
  keep_g <- apply(genes[, genes_interest, drop=FALSE],
                  2, mad, na.rm=TRUE) > 0
  keep_t <- apply(tes[, tes_interest, drop=FALSE],
                  2, mad, na.rm=TRUE) > 0
  genes_interest <- genes_interest[keep_g]
  tes_interest <- tes_interest[keep_t]
  
  stopifnot(all(length(genes_interest) > 0 & length(tes_interest) > 0))
  
  # 2a) correlate each gene to the TE‐module eigengene  
  teME     <- ME_tes[[ paste0("ME", teMod) ]]
  cor_g2T  <- bicor( genes[, genes_interest],
                     teME, use = "pairwise.complete.obs" )
  p_g2T <- corPvalueStudent(cor_g2T,
                            nSamples = length(teME))
  # cor_g2T is a vector of length = length(genes)
  gene_stats <- data.frame(
    feature = genes_interest,
    cor = as.numeric(cor_g2T),
    pval = as.numeric(p_g2T),
    qval = p.adjust(as.numeric(p_g2T), method = "BH")
  )
  gene_stats <- gene_stats[order(-abs(gene_stats$cor)), ]
  
  # 2b) correlate each TE to the gene‐module eigengene  
  geneME   <- ME_genes[[ paste0("ME", geneMod) ]]
  cor_T2g  <- bicor(tes[, tes_interest], geneME,
                    use = "pairwise.complete.obs" )
  p_T2g <- corPvalueStudent(cor_T2g,
                            nSamples = length(geneME))
  te_stats <- data.frame(
    feature = tes_interest,
    cor = as.numeric(cor_T2g),
    pval = as.numeric(p_T2g),
    qval = p.adjust(as.numeric(p_T2g), method = "BH")
  )
  te_stats <- te_stats[order(-abs(te_stats$cor)), ]
  
  # 3) return the top n of each
  list(
    topGenes = head(gene_stats, n),
    topTEs   = head(te_stats,   n)
  )
}

# ciliated
top_ciliated <- getTopPairs(geneMod = "turquoise",
                            teMod = "turquoise",
                            n = 50)
# secretory
top_secretory <- getTopPairs(geneMod = "yellow",
                             teMod = "brown",
                             n = 50)

# undiff secretory
top_sub_secretory <- getTopPairs(geneMod = "green",
                                 teMod = "yellow",
                                 n = 50)


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
  message("Significant TEs: ", nrow(x.sig))
  # plot(hist(x$summary.AUC, breaks = 40))
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
  return(x[x$FDR < 0.05 & x$summary.AUC >= 0.85,])
})

# filter to intergenic TEs here too
te_annot <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/intergenic_intronic_tes.txt.gz")
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

# remove the immune pop
immune <- colnames(te_cpm_filt.sce)[te_cpm_filt.sce$cluster %in% "6"]
to_plot <- to_plot[!rownames(to_plot) %in% immune,]

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

# heatmap of TE specific expression
library(ComplexHeatmap)

te_sig.intergenic <- unique(unlist(lapply(te_markers.sig.intergenic, function(x) rownames(x))))
sig_fte_te_exp <- logcounts(te_cpm_filt.sce)[rownames(te_cpm_filt.sce) %in% te_sig.intergenic,]

te_cpm_filt.sce.sig <- te_cpm_filt.sce[rownames(te_cpm_filt.sce) %in% te_sig.intergenic,]
table(rowData(te_cpm_filt.sce.sig)$repClass)
# LINE  LTR SINE 
# 116   99   34

# after restricting to stranded quants
# LINE  LTR SINE 
# 74   56   23 

# LINE       LTR      SINE 
# 0.4658635 0.3975904 0.1365462 

# after restricting to stranded quants
# LINE       LTR      SINE 
# 0.4836601 0.3660131 0.1503268 
# percentages are roughly the same

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

ht <- Heatmap(as.matrix(sig_fte_te_exp),
        column_split = noimmune_clust,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params) 
draw(ht)
# recover clustering of TEs for SS2 comparison
row_dend <- row_dend(ht)
row_order <- order.dendrogram(row_dend)
te_ord <- rownames(as.matrix(sig_fte_te_exp))[row_order]
write.table(te_ord,
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1_te_heatmap_order_names.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            sep = '\t')

# write out the per-cluster marker TEs from patient 1
pat1_path <- "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/"
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
            file = "~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/pat1_sig_te_markers_across_clusters.txt",
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
saveRDS(te_cpm_filt2.sce, file = "../pat2_hiseq_fte_te_cpm_sce_ens101_raw.rds")

# drag in the sce
te_cpm_filt2.sce <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/pat2/pat2_hiseq_fte_te_cpm_sce_ens101_raw.rds")

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
pat2_sce <- readRDS("~/Documents/manuscripts/storm_seq/fte_analysis/cryo_fte_hiseq_grch38_salmon_norm_se_with_embeddings.rds")
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
saveRDS(te_cpm_filt2.sce.sig, "pat2_hiseq_fte_te_cpm_sce_ens101_pat1_te_markers_with_embeddings.rds")

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
fimo <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/fimo.tsv")

# filter to significant hits
fimo.sig <- fimo[fimo$q.value < 0.05,]

# looks like zinc-finger proteins dominate 
library(dplyr)
motif_counts <- fimo.sig %>%
  filter(!is.na(motif_alt_id)) %>%           # drop rows where motif_alt_id is NA
  count(motif_alt_id, name = "abundance") %>% # tally each motif
  filter(!is.na(abundance)) %>%               # (abundance from count shouldn't be NA, but just in case)
  arrange(abundance) %>%                      # sort by count
  mutate(motif_alt_id = factor(motif_alt_id, levels = motif_alt_id)) 

# 3) Build the lollipop plot
ggplot(motif_counts, aes(x = abundance, y = motif_alt_id)) +
  # draw the “stick”
  geom_segment(aes(x = 0, xend = abundance, y = motif_alt_id, yend = motif_alt_id),
               color = "grey70") +
  # draw the “candy”
  geom_point(size = 3, color = "steelblue") +
  ylab("Enriched motif (q < 0.05)") +
  xlab("Abundance") +
  # clean theme
  theme_bw(12) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16)
  )

# make the motif
library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)

# get znf382
fimo.sig.znf382 <- fimo.sig[fimo.sig$motif_alt_id %in% "ZNF382",]

ggplot() + 
  geom_logo(seq_type = "DNA", 
            method = "prob", 
            col_scheme = "nucleotide", fimo.sig.znf382$matched_sequence) + 
  theme_logo(base_size = 24) +
  ggtitle("ZNF382 base diversity (frequency)") +
  xlab("Position") +
  scale_x_continuous(breaks = seq(1:24), labels = seq(1:24)) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )


# get znf382
fimo.sig.znf384 <- fimo.sig[fimo.sig$motif_alt_id %in% "ZNF384",]

ggplot() + 
  geom_logo(seq_type = "DNA", 
            method = "prob", 
            col_scheme = "nucleotide", fimo.sig.znf384$matched_sequence) + 
  theme_logo(base_size = 24) +
  ggtitle("ZNF384 base diversity (frequency)") +
  xlab("Position") +
  scale_x_continuous(breaks = seq(1:12), labels = seq(1:12)) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

fimo.sig.zeb1 <- fimo.sig[fimo.sig$motif_alt_id %in% "ZEB1",]

ggplot() + 
  geom_logo(seq_type = "DNA", 
            method = "prob", 
            col_scheme = "nucleotide", fimo.sig.zeb1$matched_sequence) + 
  theme_logo(base_size = 24) +
  ggtitle("ZNF384 base diversity (frequency)") +
  xlab("Position") +
  scale_x_continuous(breaks = seq(1:12), labels = seq(1:12)) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

## perform plotting of enriched GO terms

raw <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1/panther_turquoise_gene_analysis_bio.txt",
                  skip = 6)

# If the first column name is long (e.g., "GO biological process complete ..."),
# standardize column names to something usable:
if (ncol(raw) >= 8) {
  colnames(raw)[1:8] <- c("term", "ref_count", "in_list", "expected",
                          "over_under", "fold_enrichment", "p_raw", "fdr")
}

# 2) Keep only real GO rows and coerce numeric columns
df <- raw %>%
  filter(str_detect(term, "\\(GO:\\d+\\)")) %>%             # keep GO rows
  mutate(
    # Strip "<" and coerce numerics safely
    ref_count       = readr::parse_number(as.character(ref_count)),
    in_list         = readr::parse_number(as.character(in_list)),
    expected        = readr::parse_number(str_replace(as.character(expected), "<\\s*", "")),
    fold_enrichment = readr::parse_number(as.character(fold_enrichment)),
    p_raw           = readr::parse_number(as.character(p_raw)),
    fdr             = readr::parse_number(as.character(fdr)),
    over_under      = as.character(over_under),
    # Split term into name + GO ID, and make a tidy label
    go_id           = str_extract(term, "GO:\\d+"),
    go_name         = str_trim(str_remove(term, "\\s*\\(GO:\\d+\\)$")),
    label           = paste0(go_name, " (", go_id, ")")
  )

# 3) Pick the “top” enriched terms
#    Here: over-represented (+), sorted by FDR, take top N.
topN <- 10
top_over <- df %>%
  filter(over_under == "+") %>%
  arrange(fdr) %>%
  slice_head(n = topN) %>%
  # wrap long labels so they fit nicely
  mutate(label_wrapped = stringr::str_wrap(label, width = 45),
         neglogFDR = -log10(fdr))

# 4) Plot: dot plot of −log10(FDR); size = fold enrichment
ggplot(top_over,
       aes(x = fct_reorder(label_wrapped, neglogFDR),
           y = neglogFDR,
           size = fold_enrichment)) +
  geom_point(alpha = 0.85) +
  coord_flip() +
  scale_size_area(max_size = 9, name = "Fold\nenrichment") +
  labs(x = NULL,
       y = expression(-log[10]~FDR),
       title = "Top enriched GO Biological Process terms",
       subtitle = paste0("Over-represented (+); top ", topN, " by FDR")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

# --- load SEA output ---
# File must have columns: RANK, DB, ID, ALT_ID, CONSENSUS, ENR_RATIO, QVALUE, ...
sea <- read_tsv("~/Downloads/sea.tsv", show_col_types = FALSE)

# Prep stats
sea2 <- sea %>%
  mutate(
    TF   = ifelse(!is.na(ALT_ID) & ALT_ID != "", ALT_ID, ID),
    q    = suppressWarnings(as.numeric(QVALUE)),
    enr  = suppressWarnings(as.numeric(ENR_RATIO)),
    ml10 = -log10(pmax(q, 1e-300))  # stability for tiny q
  )

# --- group mapper (expanded) ---
map_group <- function(tf){
  case_when(
    # KRAB / C2H2 zinc-fingers
    str_detect(tf, "^ZNF|^ZFP|ZBTB|PATZ1") ~ "KRAB/C2H2 zinc-fingers",
    
    # NOTCH effectors (bHLH repressors)
    str_detect(tf, "^(HES1|HES2|HES5|HES6|HEY1|HEY2)$") ~ "NOTCH effectors (HES/HEY)",
    
    # Architectural / insulators
    str_detect(tf, "^(CTCF|YY1)$") ~ "Architectural factors (CTCF/YY1)",
    
    # MYC/MAX/MNT etc. (E-box bHLH)
    str_detect(tf, "^(MYC|MAX|MYCN|MNT|TFAP4|MAX::MYC)$") ~ "MYC/MAX/MXD (E-box bHLH)",
    
    # Circadian bHLH / PAR-bZIP
    str_detect(tf, "^(CLOCK|ARNTL|NPAS2|BHLHE22|BHLHE23|DBP|TEF)$") ~ "Circadian bHLH / PAR-bZIP",
    
    # HIF axis
    str_detect(tf, "^(HIF1A|EPAS1|ARNT|ARNT2|ARNT::HIF1A)$") ~ "HIF/ARNT (hypoxia)",
    
    # KLF / SP / MAZ
    str_detect(tf, "^(KLF\\d+|KLF1|KLF4|KLF5|KLF7|KLF9|KLF10|KLF12|KLF14|KLF16|SP1|SP2|SP3|SP4|SP8|MAZ)$") ~ "KLF/SP family",
    
    # FOX family (incl. FOXP/FOXA/FOXE/FOXS/FOXL)
    str_detect(tf, "^(FOX[A-Z]|FOXP\\d|FOXA\\d|FOXE\\d|FOXS\\d|FOXL\\d)$") ~ "FOX family",
    
    # TEAD
    str_detect(tf, "^(TEAD1|TEAD2|TEAD3|TEAD4)$") ~ "TEAD / YAP–TEAD",
    
    # HOX / CDX / ONECUT / OLIG / BARHL etc. (developmental homeobox)
    str_detect(tf, "^(HOX|HOXA10|HOXD9|CDX\\d|ONECUT\\d|OLIG\\d|BARHL\\d|ATOH7|MSGN1)$") ~ "Homeobox (HOX/CDX/ONECUT/OLIG)",
    
    # POU family
    str_detect(tf, "^(POU1F1|POU2F1|POU2F2|POU2F3|POU3F1|POU3F3|POU3F4|POU5F1|POU5F1B)$") ~ "POU family",
    
    # STATs
    str_detect(tf, "^(STAT\\d|Stat\\d|STAT5A|STAT5B)$") ~ "STAT family",
    
    # Nuclear receptors (RAR/RXR/PPARG/NR3C1/NR2C2/THRB/RXRG)
    str_detect(tf, "^(PPARG|Pparg::Rxra|RARB|RARG|RXRG|NR3C1|NR2C2|THRB)$") ~ "Nuclear receptors",
    
    # CREB/ATF-like bZIP (Creb3l2/Creb3l4, TFE3, TFEB, MITF)
    str_detect(tf, "^(CREB|Creb3l2|CREB3L4|TFE3|TFEB|MITF)$") ~ "bZIP (CREB/TFEB/MITF)",
    
    # ETS / ETV
    str_detect(tf, "^(ETV|ELK|ERG|ETV6)$") ~ "ETS family",
    
    # PRDM9 (meiotic hotspot setter)
    # str_detect(tf, "^PRDM9$") ~ "PRDM9 (meiotic hotspot)",
    
    TRUE ~ "Other"
  )
}

sea2 <- sea2 %>% mutate(group = map_group(TF))

# Optional: keep only reasonably enriched motifs (e.g., q <= 0.05)
sea_f <- sea2 %>% filter(is.finite(q), q <= 0.05)

# --- summarize per group ---
grp_sum <- sea_f %>%
  group_by(group) %>%
  summarize(
    best_mlog10q = max(ml10, na.rm = TRUE),     # strongest signal in group
    med_enr      = median(enr, na.rm = TRUE),   # typical fold-enrichment
    n_motifs     = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(best_mlog10q)) %>%
  mutate(group = fct_inorder(group))

# --- lollipop plot (no per-motif labels) ---
ggplot(grp_sum, aes(x = best_mlog10q, y = group)) +
  geom_segment(aes(x = 0, xend = best_mlog10q, y = group, yend = group),
               linewidth = 0.9, alpha = 0.5) +
  geom_point(aes(size = med_enr), alpha = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.10))) +
  scale_size(name = "Median fold enrichment", range = c(2.2, 6)) +
  labs(
    x = expression(paste("Best per-group  -log"[10], "(Q)")),
    y = NULL,
    title = "Motif enrichment by TF family/group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
# need to annotate ss2 TE data with cell types... fuck
ss2_all <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/new_ss2_ens102_geneannot_with_embeddings_with_celltype_annotsfte_sce.rds")

ss2_all.cell_types <- do.call(rbind, lapply(1:length(ss2_all), function(x) {
  return(data.frame(cells = colnames(ss2_all[[x]]),
                    cell_types = ss2_all[[x]]$cluster_consensus))
}))

# plot out the ss2 heatmap
ss2_tes <- readRDS("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/ss2_cpm_filt_ens101_gene_tequant_sce_sig_te_storm.rds")

# realize the matrices for plotting
counts(ss2_tes) <- as.matrix(counts(ss2_tes))
logcounts(ss2_tes) <- as.matrix(logcounts(ss2_tes))

# what is the detection looking like
min(logcounts(ss2_tes))
# 0

max(logcounts(ss2_tes))
# 11.0807

# need to go through and parse apart which cell types are which
# also label by donor

# patient33572
patient33572 <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/patient33572_file_ids.txt",
                           header = FALSE)
patient33572$patient <- "patient33572"

# patient33778
patient33778 <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/patient33778_file_ids.txt",
                           header = FALSE)
patient33778$patient <- "patient33778"

# patient34350
patient34350 <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/patient34350_file_ids.txt",
                           header = FALSE)
patient34350$patient <- "patient34350"

# patient34659
patient34659 <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/patient34659_file_ids.txt",
                           header = FALSE)
patient34659$patient <- "patient34659"

ss2_tes$patient <- ifelse(colnames(ss2_tes) %in% patient33572$V1,
                          "patient33572",
                          ifelse(colnames(ss2_tes) %in% patient33778$V1,
                                 "patient33778",
                                 ifelse(colnames(ss2_tes) %in% patient34350$V1,
                                        "patient34350",
                                        "patient34659")))

# curate clusters
ss2_tes <- ss2_tes[,colnames(ss2_tes) %in% ss2_all.cell_types$cells]
ss2_all.cell_types <- ss2_all.cell_types[ss2_all.cell_types$cells %in% colnames(ss2_tes),]
ss2_tes.m <- match(ss2_all.cell_types$cells,
                   colnames(ss2_tes))
ss2_tes <- ss2_tes[,ss2_tes.m]
all(colnames(ss2_tes) == ss2_all.cell_types$cells)
# TRUE

ss2_tes$cell_type <- ss2_all.cell_types$cell_types
ss2_tes$cell_type <- gsub("sec", "Secretory",
                          ss2_tes$cell_type)
ss2_tes$cell_type <- gsub("cil", "Ciliated",
                          ss2_tes$cell_type)
ss2_tes$cell_type <- gsub("non-epithelial", "Non-epithelial",
                          ss2_tes$cell_type)

# heatmap time
library(ComplexHeatmap)
te_ord <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/fte/pat1_te_heatmap_order_names.txt",
                     header = FALSE)

ss2_tes.m <- match(te_ord$V1,
                   rownames(ss2_tes))
ss2_tes_ord <- ss2_tes[ss2_tes.m,]
all(rownames(ss2_tes_ord) == te_ord$V1)
# TRUE

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

heatmap.anno <- HeatmapAnnotation(Cell_Type = ss2_tes_ord$cell_type,
                                  col = list(Cell_Type = c("Secretory" = "black",
                                                           "Ciliated" = "darkred",
                                                           "Non-epithelial" = "gray")),
                                     annotation_legend_param = list(
                                       Cell_Type = list(
                                         title = "Cell Type",
                                         at = unique(ss2_tes_ord$cell_type),
                                         labels = unique(ss2_tes_ord$cell_type))))

Heatmap(logcounts(ss2_tes_ord),
        column_split = ss2_tes_ord$patient,
        top_annotation = heatmap.anno,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(CPM+1)",
        heatmap_legend_param = heatmap_legend_params) 
save.image(file = "~/Documents/manuscripts/storm_seq/te_quant/fte/ss2_fte/heatmap_image_20250923.rda")
