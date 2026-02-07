## process STORM, VASA, and SSTOTAL data for gene length detection bias

## use the results out of kb python for now
## all use Ensembl 101 except for pre-processed non-UMI data which uses Ensembl 84

## libraries
library(rtracklayer)
library(GenomicFeatures)
library(eisaR)
library(scuttle)
library(scran)
library(Matrix)

## pull in the gtf and extract and summarize exons
ens101.gtf <- rtracklayer::import("~/Documents/kb_python_storm_test/ens101_ref/Homo_sapiens.GRCh38.101.ercc92patched.gtf.gz")
ens84.gtf <- rtracklayer::import("~/Documents/manuscripts/mini_scrna_2020/gene_length_bias/GSE63818-GPL16791/ref/Homo_sapiens.GRCh38.84.gtf")

## filter to standard chromosomes (exclude MT and ERCCs)
ens101.gtf.stdchroms <- keepSeqlevels(ens101.gtf,
                                      c(as.character(1:22), "X", "Y"),
                                      pruning.mode = "coarse")
ens84.gtf.stdchroms <- keepSeqlevels(ens84.gtf,
                                      c(as.character(1:22), "X", "Y"),
                                      pruning.mode = "coarse")

## remove ribo genes
ens101.gtf.stdchroms.noribo <- ens101.gtf.stdchroms[!ens101.gtf.stdchroms$gene_biotype %in% 
                                                      c("rRNA", "rRNA_pseudogene",
                                                        "Mt_rRNA"),]
ens84.gtf.stdchroms.noribo <- ens84.gtf.stdchroms[!ens84.gtf.stdchroms$gene_biotype %in% 
                                                      c("rRNA", "rRNA_pseudogene",
                                                        "Mt_rRNA"),]

## make the TxDb
ens101.txdb <- makeTxDbFromGRanges(ens101.gtf.stdchroms.noribo)
ens84.txdb <- makeTxDbFromGRanges(ens84.gtf.stdchroms.noribo)

# get exons and summarize to gene lengths for binning
ens101.exons <- getRegionsFromTxDb(ens101.txdb,
                                   strandedData = TRUE)
ens84.exons <- getRegionsFromTxDb(ens84.txdb,
                                   strandedData = TRUE)
saveRDS(ens101.exons,
        file = "~/Documents/kb_python_storm_test/ens101_ref/eisaR_ens101_exons.rds")
saveRDS(ens84.exons,
        file = "~/Documents/manuscripts/mini_scrna_2020/gene_length_bias/GSE63818-GPL16791/ref/eisaR_ens84_exons.rds")
ens101.genelengths <- sort(ens101.exons$genebodies)
ens101.genelengths.flat <- width(ens101.genelengths)
names(ens101.genelengths.flat) <- names(ens101.genelengths)

ens84.genelengths <- sort(ens84.exons$genebodies)
ens84.genelengths.flat <- width(ens84.genelengths)
names(ens84.genelengths.flat) <- names(ens84.genelengths)

## modified function from https://github.com/Oshlack/GeneLengthBias-scRNASeq/blob/master/Klein-hK562-2015.Rmd
binnedGenes <- function(gr) {
  q<-quantile(sqrt(width(ens101.genelengths)),probs=seq(0.1,1,0.1))
  decile <- rep(NA,length(gr))
  decile[sqrt(width(ens101.genelengths))<=q[1]] <- 1
  for(i in 2:10) decile[sqrt(width(gr))>q[i-1] & sqrt(width(gr))<=q[i]] <- i
  return(decile)
} 

# generate deciles of gene lengths by exons
gene_bins <- binnedGenes(ens101.genelengths)

# make the names of each bin the gene name
names(gene_bins) <- names(ens101.genelengths)

# save the bins
saveRDS(gene_bins,
        file = "~/Documents/kb_python_storm_test/ens101_ref/decile_ens101_gene_length_bins.rds")

# fix the gene rownames
fix_rownames <- function(txis, stub="^ENS", sep="\\.", idx=1) { 
  
  fixable <- grep(stub, rownames(txis))
  if (length(fixable) < 1) { 
    message("No fixable rownames. Returning unaltered.")
    return(txis)
  }
  rownames(txis)[fixable] <- 
    sapply(strsplit(rownames(txis)[fixable], sep), `[`, idx)
  message("Fixed ", length(fixable), " rownames.") 
  return(txis) 
  
}

## preprocess and QC data out of kb
read_kb_gene <- function(res_dir, mtx_pattern = "cells_x_genes.mtx",
                         genes_file = "genes.txt",
                         cell_name=NULL) {
  ## list the files and find the mtx
  mtx <- list.files(res_dir,
                    pattern = mtx_pattern,
                    full.names = TRUE)
  # read in
  mat.mm <- readMM(mtx)
  # transpose
  mat.mm <- Matrix::t(mat.mm)
  mat.mm <- as(mat.mm, "dgCMatrix")
  genes <- list.files(res_dir,
                      pattern = genes_file,
                      full.names = TRUE)
  # read in
  genes.tab <- read.delim(genes, header = FALSE)
  if (is.null(cell_name)) {
    colnames(mat.mm) <- "cell"
  } else {
    colnames(mat.mm) <- cell_name
  }
  rownames(mat.mm) <- genes.tab$V1
  return(mat.mm)
}

gatherMMFiles <- function(path, file = "cells_x_genes.mtx",
                          dir_suffix = "_kb_quant/counts_unfiltered") {
  mmf <- list.files(path = path,
                    pattern = file,
                    recursive = TRUE,
                    full.names = TRUE)
  mmf <- dirname(mmf)
  fgsub <- gsub(paste0(path, "/"),
                "",
                mmf)
  names(mmf) <- gsub(dir_suffix,
                     "",
                     fgsub)
  return(mmf)
}

ingestKB <- function(mm_files, file = "cells_x_genes.mtx",
                     genes_file = "genes.txt") {
  stopifnot(!is.null(names(mm_files)))
  raw_counts <- lapply(names(mm_files), function(x) {
    kb_mm <- read_kb_gene(mm_files[x], cell_name = x,
                          mtx_pattern = file,genes_file = genes_file)
  })
  raw_counts.mat <- do.call(cbind, raw_counts)
  return(SingleCellExperiment(assays = list(counts = raw_counts.mat)))
}

# need to specify pos and neg control wells in final map
# need to specify a named vector of cell types
remapStorm <- function(sce, remap_file = NULL, header = FALSE,
                       pos_control_wells = NULL,
                       neg_control_wells = NULL,
                       cell_types = NULL) {
  if (is.null(remap_file)) stop("Need to specify a two-column well remapping file.")
  storm_remap <- read.delim(remap_file,
                            header = header)
  # match
  message("Remapping wells.")
  storm_remap.match <- match(colnames(sce),
                             storm_remap$V1)
  storm_remap <- storm_remap[storm_remap.match,]
  
  #careful
  stopifnot(all(storm_remap$V1 == colnames(sce)))
  # TRUE
  
  # re-map the names/wells
  colnames(sce) <- storm_remap$V2
  
  if (!is.null(pos_control_wells)) {
    message("Filtering out positive controls.")
    sce <- sce[,!colnames(sce) %in% pos_control_wells]
  }
  
  if (!is.null(neg_control_wells)) {
    message("Filtering out negative controls.")
    sce <- sce[,!colnames(sce) %in% neg_control_wells]
  }
  
  if (!is.null(cell_types)) {
    message("Annotating cell types.")
    cell_types.match <- match(colnames(sce),
                              cell_types)
    cell_types <- cell_types[cell_types.match]
    stopifnot(all(cell_types == colnames(sce)))
    sce$cell_type <- names(cell_types)
    message("Cell types:")
    print(table(sce$cell_type))
  }
  return(sce)
}

splitSpikeIns <- function(sce, spike_in = "ERCC-",
                          altExp_slot_name = "ERCC") {
  # pull out the spike ins
  spike.sce <- sce[grep(spike_in, rownames(sce)),]
  if (nrow(spike.sce) == 0) {
    message("No spike-ins found. Returning object as is.")
    return(sce)
  }
  message("Splitting out the spike-ins.")
  altExp(sce, altExp_slot_name) <- spike.sce
  sce <- sce[grep(spike_in, rownames(sce), invert = TRUE),]
  
  # re-calculate num genes
  message("Recalculating number of genes expressed.")
  sce$NumGenesExpressed <- colSums2(counts(sce) > 0)
  return(sce)
}

quickQC <- function(sce, sub.fields = NULL) {
  sce <- addPerCellQCMetrics(sce)
  filt_reasons <- perCellQCFilters(sce,
                                   sub.fields = sub.fields)
  message("Filter reasons:")
  print(colSums(as.matrix(filt_reasons)))
  
  # filter
  message("Filtering poor quality cells.")
  sce.filt <- sce[,!filt_reasons$discard]
  return(sce.filt)
}

readDepthFilterCells <- function(sce,
                                 reads_file = "read_counts.txt",
                                 min_threshold = 1e5,
                                 file_suffix = ".100k.fq.gz") {
  # read in file
  read_depth <- read.delim(reads_file,
                           header = FALSE)
  # find cells with sufficient depth
  read_depth.filt <- read_depth[read_depth[,2] >= (min_threshold*4),]
  # filter
  read_depth.filt[,1] <- gsub(file_suffix, "", read_depth.filt[,1])
  sce.filt <- sce[,colnames(sce) %in% read_depth.filt[,1]]
  return(sce.filt)
}


## SSTOTAL
# HEK293T
sstotal_mm_files <- gatherMMFiles("/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/sstotal",
                                  file = "matrix.abundance.gene.mtx",
                                  dir_suffix = "_kb_tcc/quant_unfiltered")

sstotal_sce <- ingestKB(sstotal_mm_files, file = "matrix.abundance.gene.mtx")
sstotal_sce <- splitSpikeIns(sstotal_sce)
sstotal_sce <- fix_rownames(sstotal_sce)

# filter to 100k cells
sstotal_sce <- readDepthFilterCells(sstotal_sce,
                                    reads_file = "/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/sstotal/read_counts.txt",
                                    file_suffix = "_trimmed.100k.fq.gz")
sstotal_sce <- quickQC(sstotal_sce)

# filter genes
gene_filter <- (rowSums(counts(sstotal_sce)==0)/ncol(sstotal_sce)) <= 0.9
sstotal_sce <- sstotal_sce[gene_filter,]

# normalize
sstotal_sce <- logNormCounts(sstotal_sce)

# subset
sstotal_sce <- sstotal_sce[rownames(sstotal_sce) %in% names(ens101.genelengths.flat),]

# calculate FPKM using above gene lengths
sstotal_gl_match <- match(rownames(sstotal_sce),
                          names(ens101.genelengths.flat))

sstotal_gl <- ens101.genelengths.flat[sstotal_gl_match]

# check
all(names(sstotal_gl) == rownames(sstotal_sce))
# TRUE

rowData(sstotal_sce)$gene_length <- sstotal_gl

assay(sstotal_sce, "log_fpkm") <- log2(calculateFPKM(sstotal_sce,
                                                      lengths = rowData(sstotal_sce)$gene_length) + 1)

## add on the bin number
propZ_genes <- rowSums(counts(sstotal_sce)==0)/ncol(sstotal_sce)
q<-quantile(sqrt(rowData(sstotal_sce)$gene_length),probs=seq(0.1,1,0.1))
decile <- rep(NA,nrow(sstotal_sce))
decile[sqrt(rowData(sstotal_sce)$gene_length)<=q[1]] <- 1
for(i in 2:10) decile[sqrt(rowData(sstotal_sce)$gene_length)>q[i-1] & sqrt(rowData(sstotal_sce)$gene_length)<=q[i]] <- i


# plot data
sstotal_gene_length_bias_data <- data.frame(log_counts = rowMeans(logcounts(sstotal_sce)),
                                            log_fpkm = rowMeans(assay(sstotal_sce, "log_fpkm")),
                                            prop_zero = propZ_genes,
                                            bins = factor(decile))

library(ggplot2)

ggplot(sstotal_gene_length_bias_data,
       aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,8))

ggplot(sstotal_gene_length_bias_data,
       aes(x = bins, y = prop_zero, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(0,1)

ggplot(sstotal_gene_length_bias_data,
       aes(x = bins, y = log_fpkm, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,14))

## VASA
# HEK293T
vasa_mm_files <- gatherMMFiles("/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/vasa",
                                  file = "matrix.abundance.gene.mtx",
                                  dir_suffix = "_kb_tcc/quant_tcc")

vasa_sce <- ingestKB(vasa_mm_files, file = "matrix.abundance.gene.mtx")
vasa_sce <- splitSpikeIns(vasa_sce)
vasa_sce <- fix_rownames(vasa_sce)

# split into HEK cells
colnames(vasa_sce)[1:768] <- c(paste0(colnames(vasa_sce)[1:192], "_K562"),
                               paste0(colnames(vasa_sce)[193:384], "_HEK293T"),
                               paste0(colnames(vasa_sce)[385:576], "_K562"),
                               paste0(colnames(vasa_sce)[577:768], "_HEK293T"))

p1_neg_controls <- c("VAI-MR-v001_HNVCGBGXN_S3_009_K562", "VAI-MR-v001_HNVCGBGXN_S3_047_K562",
                     "VAI-MR-v001_HNVCGBGXN_S3_108_K562", "VAI-MR-v001_HNVCGBGXN_S3_123_K562",
                     "VAI-MR-v001_HNVCGBGXN_S3_207_HEK293T", "VAI-MR-v001_HNVCGBGXN_S3_318_HEK293T",
                     "VAI-MR-v001_HNVCGBGXN_S3_334_HEK293T", "VAI-MR-v001_HNVCGBGXN_S3_379_HEK293T")
p2_neg_controls <- c("VAI-MR-v002_HNVCGBGXN_S4_009_K562", "VAI-MR-v002_HNVCGBGXN_S4_047_K562",
                     "VAI-MR-v002_HNVCGBGXN_S4_108_K562", "VAI-MR-v002_HNVCGBGXN_S4_123_K562",
                     "VAI-MR-v002_HNVCGBGXN_S4_207_HEK293T", "VAI-MR-v002_HNVCGBGXN_S4_318_HEK293T",
                     "VAI-MR-v002_HNVCGBGXN_S4_334_HEK293T", "VAI-MR-v002_HNVCGBGXN_S4_379_HEK293T")
p1_neg_controls %in% colnames(vasa_sce)
# TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
p2_neg_controls %in% colnames(vasa_sce)
# TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#filter
vasa_sce <- vasa_sce[,!colnames(vasa_sce) %in% c(p1_neg_controls,
                                                 p2_neg_controls)]
vasa_sce <- vasa_sce[,grep("HEK293T", colnames(vasa_sce))]

# filter to 100k read depth cells
# need to fix the colnames to remove _HEK293T
colnames(vasa_sce) <- gsub("_HEK293T",
                           "",
                           colnames(vasa_sce)) 
vasa_sce <- readDepthFilterCells(vasa_sce,
                                 reads_file = "/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/vasa/read_counts.txt",
                                 file_suffix = "_cbc_trimmed_homoATCG.fixed.100k.fq.gz")

vasa_sce <- quickQC(vasa_sce)

# filter genes
gene_filter <- (rowSums(counts(vasa_sce)==0)/ncol(vasa_sce)) <= 0.9
vasa_sce <- vasa_sce[gene_filter,]

# normalize
vasa_sce <- logNormCounts(vasa_sce)

# subset
vasa_sce <- vasa_sce[rownames(vasa_sce) %in% names(ens101.genelengths.flat),]

# calculate FPKM using above gene lengths
vasa_gl_match <- match(rownames(vasa_sce),
                          names(ens101.genelengths.flat))

vasa_gl <- ens101.genelengths.flat[vasa_gl_match]

# check
all(names(vasa_gl) == rownames(vasa_sce))
# TRUE

rowData(vasa_sce)$gene_length <- vasa_gl

assay(vasa_sce, "log_fpkm") <- log2(calculateFPKM(vasa_sce,
                                                     lengths = rowData(vasa_sce)$gene_length) + 1)

## add on the bin number
propZ_genes <- rowSums(counts(vasa_sce)==0)/ncol(vasa_sce)
q<-quantile(sqrt(rowData(vasa_sce)$gene_length),probs=seq(0.1,1,0.1))
decile <- rep(NA,nrow(vasa_sce))
decile[sqrt(rowData(vasa_sce)$gene_length)<=q[1]] <- 1
for(i in 2:10) decile[sqrt(rowData(vasa_sce)$gene_length)>q[i-1] & sqrt(rowData(vasa_sce)$gene_length)<=q[i]] <- i

# plot data
vasa_gene_length_bias_data <- data.frame(log_counts = rowMeans(logcounts(vasa_sce)),
                                            log_fpkm = rowMeans(assay(vasa_sce, "log_fpkm")),
                                            prop_zero = propZ_genes,
                                            bins = factor(decile))


ggplot(vasa_gene_length_bias_data,
       aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d()

ggplot(vasa_gene_length_bias_data,
       aes(x = bins, y = prop_zero, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(0,1)

ggplot(vasa_gene_length_bias_data,
       aes(x = bins, y = log_fpkm, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d()

## STORM
# HEK293T
storm_mm_files <- gatherMMFiles("/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/storm",
                               file = "matrix.abundance.gene.mtx",
                               dir_suffix = "_kb_tcc/quant_tcc")

storm_sce <- ingestKB(storm_mm_files, file = "matrix.abundance.gene.mtx")
storm_sce <- splitSpikeIns(storm_sce)
storm_sce <- fix_rownames(storm_sce)
pos_control <- c("H2",
                 "G7")
neg_control <- c("L2",
                 "I7")
k562_cells <- c(paste0("A", 1:24),
                paste0("D", 1:24),
                paste0("F", 1:24),
                paste0("K", 1:24),
                paste0("P", 1:24))
names(k562_cells) <- rep("K562", length(k562_cells))

hek293t_cells <- c(paste0("B", 1:24),
                   paste0("E", 1:24),
                   paste0("G", 1:6),
                   paste0("G", 8:24),
                   paste0("I", 1:6),
                   paste0("I", 8:24),
                   paste0("L", c(1,3:24)),
                   paste0("N", 1:24))
names(hek293t_cells) <- rep("HEK293T", length(hek293t_cells))

rmg2_cells <- c(paste0("C", 1:24),
                paste0("H", c(1,3:24)),
                paste0("J", 1:24),
                paste0("M", 1:24),
                paste0("O", 1:24))
names(rmg2_cells) <- rep("RMG2", length(rmg2_cells))
all_celltypes <- c(k562_cells,
                   hek293t_cells,
                   rmg2_cells)
storm_sce <- remapStorm(storm_sce, remap_file = "~/Documents/kb_python_storm_test/well_map.txt",
                        pos_control_wells = pos_control, neg_control_wells = neg_control,
                        cell_types = all_celltypes)
storm_sce <- storm_sce[,storm_sce$cell_type %in% "HEK293T"]

# filter to cells with 100k reads
# this will be tricky since they need to be remapped...
remap_file = "~/Documents/kb_python_storm_test/well_map.txt"
read_depth = read.delim("/Volumes/projects/shen/projects/SHEH_20220829_scRNA/analysis/kb_python/storm/read_counts.txt",
                        header = FALSE)

# helper
rmAfterNthDelim <- function(vec, delim="-", n=3) {
  return(substr(vec, 1, sapply(gregexpr(delim, vec), "[", n) - 1))
}

read_depth$V1 <- rmAfterNthDelim(read_depth$V1, delim = "_", n = 1)
read_depth <- read_depth[!duplicated(read_depth$V1),]
read_depth.filt <- read_depth[read_depth$V2 >= 4e5,]
storm_sce <- storm_sce[,colnames(storm_sce) %in% read_depth.filt$V1]

storm_sce <- quickQC(storm_sce, sub.fields = "altexps_ERCC_percent")

# filter genes
gene_filter <- (rowSums(counts(storm_sce)==0)/ncol(storm_sce)) <= 0.9
storm_sce <- storm_sce[gene_filter,]

# normalize
storm_sce <- computeSpikeFactors(storm_sce, "ERCC")
storm_sce <- logNormCounts(storm_sce, size.factors = storm_sce$sizeFactor)

# subset
storm_sce <- storm_sce[rownames(storm_sce) %in% names(ens101.genelengths.flat),]

# calculate FPKM using above gene lengths
storm_gl_match <- match(rownames(storm_sce),
                       names(ens101.genelengths.flat))

storm_gl <- ens101.genelengths.flat[storm_gl_match]

# check
all(names(storm_gl) == rownames(storm_sce))
# TRUE

rowData(storm_sce)$gene_length <- storm_gl

assay(storm_sce, "log_fpkm") <- log2(calculateFPKM(storm_sce,
                                                   lengths = rowData(storm_sce)$gene_length,
                                                   size.factors = sizeFactors(storm_sce)) + 1)

## add on the bin number
propZ_genes <- rowSums(counts(storm_sce)==0)/ncol(storm_sce)
q<-quantile(sqrt(rowData(storm_sce)$gene_length),probs=seq(0.1,1,0.1))
decile <- rep(NA,nrow(storm_sce))
decile[sqrt(rowData(storm_sce)$gene_length)<=q[1]] <- 1
for(i in 2:10) decile[sqrt(rowData(storm_sce)$gene_length)>q[i-1] & sqrt(rowData(storm_sce)$gene_length)<=q[i]] <- i

# plot data
storm_gene_length_bias_data <- data.frame(log_counts = rowMeans(logcounts(storm_sce)),
                                         log_fpkm = rowMeans(assay(storm_sce, "log_fpkm")),
                                         prop_zero = propZ_genes,
                                         bins = factor(decile))


ggplot(storm_gene_length_bias_data,
       aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d()

ggplot(storm_gene_length_bias_data,
       aes(x = bins, y = prop_zero, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(0,1)

ggplot(storm_gene_length_bias_data,
       aes(x = bins, y = log_fpkm, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d()


# some comparisons with non-UMI and UMI containing to show what the trends are

## sanity check to make sure things are working properly
## It is sane and works with the gene lengths
## big note though to filter genes with >90% zeroes at the start

## UMI example
library(edgeR)
library(RColorBrewer)
library(org.Hs.eg.db)
library(limma)
counts.klein <- read.csv("~/Downloads/GSM1599500_K562_cells.csv",stringsAsFactors = FALSE,row.names=1)
colnames(counts.klein) <- paste("Cell",1:ncol(counts.klein),sep="")

# Separate ercc and endogenous genes
ercc <- grep("ERCC", rownames(counts.klein))
counts.end <- counts.klein[-ercc,]
counts.ercc <- counts.klein[ercc,]

# Calculate dropout and library size
dropout <- colSums(counts.end==0)/nrow(counts.end)
lib.size <- colSums(counts.klein)
lib.size.ercc <- colSums(counts.ercc)
lib.size.end <- colSums(counts.end)

keep1 <- dropout<0.85 & lib.size.end>10000 & lib.size.ercc/lib.size < 0.01
counts.keep <- counts.end[,keep1]
dim(counts.keep)


## important to note - genes here were filtered out that had >90% zeroes
propZ_genes <- rowSums(counts.keep==0)/ncol(counts.keep)
counts.keep <- counts.keep[propZ_genes<=0.9,]
dim(counts.keep)

# make into SingleCellExperiment obj
y <- SingleCellExperiment(assay = list(counts=counts.keep))

#annotation
symbol <- toTable(org.Hs.egSYMBOL)
m <- match(rownames(y),symbol$symbol)
ann <- data.frame(Original_ID=rownames(y),symbol[m,])
rownames(ann) <- rownames(y)
ens <- toTable(org.Hs.egENSEMBL)
m <- match(ann$gene_id,ens$gene_id)
ann$ensembl_id <- ens$ensembl_id[m]
chr <- toTable(org.Hs.egCHR)
m <- match(ann$gene_id,chr$gene_id)
ann$chr <- chr$chromosome[m]
genename <- toTable(org.Hs.egGENENAME)
m <- match(ann$gene_id,genename$gene_id)
ann$genename <- genename$gene_name[m]
# swap out for our own lengths
# m <- match(ann$ensembl_id,hg38.length$EnsID)
m <- match(ann$ensembl_id,
           names(ens101.genelengths.flat))
# ann$length <- hg38.length$Length[m]
ann$length <- ens101.genelengths.flat[m]
rowData(y) <- ann

y <- y[complete.cases(rowData(y)$length),]
mito <- grep("mitochondrial", rowData(y)$genename)
ribo <- grep("ribosomal", rowData(y)$genename)
chrm <- grep("MT", rowData(y)$chr)
junk <- unique(c(mito,ribo,chrm))
length(junk)
y <- y[-junk,]
y$lib.size <- colSums(counts(y))

y <- logNormCounts(y)
assay(y, "log_fpkm") <- log2(calculateFPKM(y,
                                           lengths = rowData(y)$length,
                                           size.factors = sizeFactors(y)) + 1)
propZ_genes <- rowSums(counts(y)==0)/ncol(y)
q<-quantile(sqrt(rowData(y)$length),probs=seq(0.1,1,0.1))
decile <- rep(NA,nrow(y))
decile[sqrt(rowData(y)$length)<=q[1]] <- 1
for(i in 2:10) decile[sqrt(rowData(y)$length)>q[i-1] & sqrt(rowData(y)$length)<=q[i]] <- i

# plot data
tenx_gene_length_bias_data <- data.frame(log_counts = rowMeans(logcounts(y)),
                                          log_fpkm = rowMeans(assay(y, "log_fpkm")),
                                          prop_zero = propZ_genes,
                                          bins = factor(decile))


ggplot(tenx_gene_length_bias_data,
       aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,8))

ggplot(tenx_gene_length_bias_data,
       aes(x = bins, y = prop_zero, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(0,1)

ggplot(tenx_gene_length_bias_data,
       aes(x = bins, y = log_fpkm, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,10))


## Non-UMI example
library(tximport)

#make the tx2g
t2g <- read.delim("~/Documents/manuscripts/mini_scrna_2020/gene_length_bias/GSE63818-GPL16791/ref/ens_84_t2g.txt",
                  header = FALSE)

salmon_files <- list.files("~/Documents/manuscripts/mini_scrna_2020/gene_length_bias/GSE63818-GPL16791/salmon",
                           pattern = "quant.sf",
                           full.names = TRUE,
                           recursive = TRUE)
names(salmon_files) <- paste0("cell_", 1:length(salmon_files))
salmon.txi <- tximport(salmon_files,
                       type = "salmon",
                       tx2gene = t2g,
                       importer=read.delim)

counts <- salmon.txi$counts
tpm <- salmon.txi$abundance

# keep PGCs
keep.cells <- c(1:233)
counts.pgc <- counts[,keep.cells]
tpm.pgc <- tpm[,keep.cells]

par(mar=c(5,4,2,2))
par(mfrow=c(1,2))
dropout <- colSums(counts.pgc==0)/nrow(counts.pgc)
lib.size <- colSums(counts.pgc)
plot(dropout,lib.size)
abline(v=0.85,h=500000,lty=2,col=2)
points(dropout[dropout>0.85],lib.size[dropout>0.85],col=2,pch=16)
mycol <- rep(NA,ncol(counts.pgc))
mycol[dropout>0.85] <- 2
mycol[dropout<=0.85] <- 1
plotMDS(DGEList(counts.pgc),pch=16,gene.selection = "common",col=mycol)

# filter
counts.keep <- counts.pgc[,dropout<=0.85]
tpm.keep <- tpm.pgc[,dropout<=0.85]
dim(tpm.keep)

propZ_genes <- rowSums(counts.keep==0)/ncol(counts.keep)
counts.keep <- counts.keep[propZ_genes<=0.9,]
tpm.keep <- tpm.keep[propZ_genes<=0.9,]
dim(counts.keep)
dim(tpm.keep)
table(rownames(counts.keep)==rownames(tpm.keep))

y <- SingleCellExperiment(assay = list(counts=counts.keep))

y <- fix_rownames(y)

#annotation
ens <- toTable(org.Hs.egENSEMBL)
m <- match(rownames(y),ens$ensembl_id)
ann <- data.frame(Original_ID=rownames(y),ens[m,])

symbol <- toTable(org.Hs.egSYMBOL)
m <- match(ann$gene_id,symbol$gene_id)

ann <- cbind(ann, symbol[m,])
rownames(ann) <- rownames(y)

chr <- toTable(org.Hs.egCHR)
m <- match(ann$gene_id,chr$gene_id)
ann$chr <- chr$chromosome[m]
genename <- toTable(org.Hs.egGENENAME)
m <- match(ann$gene_id,genename$gene_id)
ann$genename <- genename$gene_name[m]
# swap out for our own lengths
# m <- match(ann$ensembl_id,hg38.length$EnsID)
m <- match(ann$ensembl_id,
           names(ens84.genelengths.flat))
# ann$length <- hg38.length$Length[m]
ann$length <- ens84.genelengths.flat[m]
rowData(y) <- ann

y <- y[complete.cases(rowData(y)$length),]
mito <- grep("mitochondrial", rowData(y)$genename)
ribo <- grep("ribosomal", rowData(y)$genename)
chrm <- grep("MT", rowData(y)$chr)
junk <- unique(c(mito,ribo,chrm))
length(junk)
y <- y[-junk,]
y$lib.size <- colSums(counts(y))

y <- logNormCounts(y)
assay(y, "log_fpkm") <- log2(calculateFPKM(y,
                                           lengths = rowData(y)$length,
                                           size.factors = sizeFactors(y)) + 1)
propZ_genes <- rowSums(counts(y)==0)/ncol(y)
q<-quantile(sqrt(rowData(y)$length),probs=seq(0.1,1,0.1))
decile <- rep(NA,nrow(y))
decile[sqrt(rowData(y)$length)<=q[1]] <- 1
for(i in 2:10) decile[sqrt(rowData(y)$length)>q[i-1] & sqrt(rowData(y)$length)<=q[i]] <- i

# plot data
ss_gene_length_bias_data <- data.frame(log_counts = rowMeans(logcounts(y)),
                                         log_fpkm = rowMeans(assay(y, "log_fpkm")),
                                         prop_zero = propZ_genes,
                                         bins = factor(decile))


ggplot(ss_gene_length_bias_data,
       aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,8))

ggplot(ss_gene_length_bias_data,
       aes(x = bins, y = prop_zero, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(0,1)

ggplot(ss_gene_length_bias_data,
       aes(x = bins, y = log_fpkm, group = bins, fill = bins)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  ylim(c(0,10))

# make a "mega counts object" to plot

sstotal_gene_length_bias_data$tech <- "SMART-seq Total"
median(sstotal_gene_length_bias_data$log_counts)
# 0.4666984

vasa_gene_length_bias_data$tech <- "VASA-seq"
median(vasa_gene_length_bias_data$log_counts)
# 0.7247008

storm_gene_length_bias_data$tech <- "STORM-seq"
median(storm_gene_length_bias_data$log_counts)
# 1.132218

tenx_gene_length_bias_data$tech <- "UMI-Droplet"
median(tenx_gene_length_bias_data$log_counts)
# 0.5077706

ss_gene_length_bias_data$tech <- "Non-UMI NEBNext Ultra"
median(ss_gene_length_bias_data$log_counts)
# 2.661144

all_data <- rbind(sstotal_gene_length_bias_data,
                  vasa_gene_length_bias_data,
                  storm_gene_length_bias_data,
                  tenx_gene_length_bias_data,
                  ss_gene_length_bias_data)


median_lines <- data.frame(tech = c("SMART-seq Total",
                                    "VASA-seq",
                                    "STORM-seq",
                                    "UMI-Droplet",
                                    "Non-UMI NEBNext Ultra"),
                           lines = c(0.4666984, 0.7247008,
                                   1.132218, 0.5077706,
                                   2.661144))

ggplot(all_data, aes(x = bins, y = log_counts, group = bins, fill = bins)) +
  geom_boxplot() +
  #geom_hline(data = median_lines, aes(yintercept = lines), color = "red") +
  scale_fill_manual(values = c(rep("white", 10))) +
  facet_wrap(~tech, scales = "free") +
  theme_bw(12) +
  theme(
    legend.position = "None"
  )
