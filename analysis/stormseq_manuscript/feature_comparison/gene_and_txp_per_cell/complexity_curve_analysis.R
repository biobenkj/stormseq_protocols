# complexity curves for storm, sstotal, vasa, and ss3

## Main functions

library(Matrix)
library(SingleCellExperiment)

# write a small function to annotate the matrix market format
read_kb_gene <- function(res_dir, cell_name=NULL) {
  ## list the files and find the mtx
  mtx <- list.files(res_dir,
                    pattern = "matrix.abundance.gene.mtx",
                    full.names = TRUE)
  # read in
  mat.mm <- readMM(mtx)
  # transpose
  mat.mm <- Matrix::t(mat.mm)
  mat.mm <- as(mat.mm, "dgCMatrix")
  genes <- list.files(res_dir,
                      pattern = "genes.txt",
                      full.names = TRUE)
  # read in
  genes.tab <- read.delim(genes, header = FALSE)
  if (is.null(cell_name)) {
    colnames(mat.mm) <- paste0("cell_", 1:ncol(mat.mm))
  } else {
    colnames(mat.mm) <- cell_name
  }
  rownames(mat.mm) <- genes.tab$V1
  return(mat.mm)
}

read_kb_txp <- function(res_dir, cell_name=NULL) {
  ## list the files and find the mtx
  mtx <- list.files(res_dir,
                    pattern = "matrix.abundance.mtx",
                    full.names = TRUE)
  # read in
  mat.mm <- readMM(mtx)
  # transpose
  mat.mm <- Matrix::t(mat.mm)
  mat.mm <- as(mat.mm, "dgCMatrix")
  genes <- list.files(res_dir,
                      pattern = "transcripts.txt",
                      full.names = TRUE)
  # read in
  genes.tab <- read.delim(genes, header = FALSE)
  if (is.null(cell_name)) {
    colnames(mat.mm) <- paste0("cell_", 1:ncol(mat.mm))
  } else {
    colnames(mat.mm) <- cell_name
  }
  rownames(mat.mm) <- genes.tab$V1
  return(mat.mm)
}

gather_kb_res <- function(path) {
  ## first gather the directories
  analysis_dirs <- list.dirs(path = path,
                             full.names = TRUE,
                             recursive = FALSE)
  ## now dip into each dir and import/summarize the gene counts
  kb_res_lst <- lapply(analysis_dirs, function(x) {
    message("Working on dir: ", x)
    kb_gene_count_files <- list.files(x,
                                      pattern = "matrix.abundance.gene.mtx",
                                      recursive = TRUE,
                                      full.names = TRUE)
    kb_gene_count_files <- dirname(kb_gene_count_files)
    # rename
    first_gsub <- gsub(x, "", kb_gene_count_files)
    names(kb_gene_count_files) <- gsub("_kb_tcc/quant_tcc",
                                       "", first_gsub)
    message("Extracting gene counts...")
    kb_gene_counts <- lapply(names(kb_gene_count_files), function(y) {
      kb_mm <- read_kb_gene(kb_gene_count_files[y],
                            y)
      return(kb_mm)
    })
    kb_txp_counts <- lapply(names(kb_gene_count_files), function(y) {
      kb_mm <- read_kb_txp(kb_gene_count_files[y],
                            y)
      return(kb_mm)
    })
    message("Merging gene counts...")
    kb_gene_counts <- do.call(cbind,
                              kb_gene_counts)
    kb_txp_counts <- do.call(cbind,
                              kb_txp_counts)
    sce.gene.raw <- SingleCellExperiment(assays = list(counts = kb_gene_counts))
    sce.txp.raw <- SingleCellExperiment(assays = list(counts = kb_txp_counts))
    return(list(gene = sce.gene.raw,
                txp = sce.txp.raw))
  })
  return(kb_res_lst)
}

# extract for storm
storm_kb_res <- gather_kb_res(getwd())
names(storm_kb_res) <- list.dirs(full.names = FALSE,
                                 recursive = FALSE)
saveRDS(storm_kb_res, file = "storm_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds") 
storm_kb_res <- readRDS("storm/storm_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")

# now need to go through and pull out the gene and txp counts that are not erccs
storm_kb_res_gene_counts <- lapply(storm_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  txp_counts <- x[["txp"]]
  txp_counts <- counts(txp_counts[grep("ERCC-",
                                       invert = TRUE,
                                       rownames(txp_counts)),])
  detected_genes <- colSums2(gene_counts >= 1)
  umi_frags <- colSums2(gene_counts)
  detected_txps <- colSums2(txp_counts >= 1)
  # assembl
  return(data.frame(detected_genes = detected_genes,
                    detected_txps = detected_txps,
                    umi_frags = umi_frags,
                    row.names = cells))
})
saveRDS(storm_kb_res_gene_counts, file = "storm_gene_txp_counts_ens101_with_erccs_for_plotting.rds") 
storm_kb_res_gene_counts <- readRDS("storm/storm_gene_txp_counts_ens101_with_erccs_for_plotting.rds")

# remap the wells to the proper place and remove controls
storm_remap <- read.delim("storm/well_map.txt",
                          header = FALSE)

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

storm_kb_res_gene_counts.remap <- lapply(storm_kb_res_gene_counts, function(x) {
  storm_remap.match <- match(rownames(x),
                             storm_remap$V1)
  storm_remap <- storm_remap[storm_remap.match,]
  stopifnot(all(storm_remap$V1 == rownames(x)))
  rownames(x) <- storm_remap$V2
  return(x)
})

# annotate cell types
hek293t_cells <- c(paste0("B", 1:24),
                   paste0("E", 1:24),
                   paste0("G", 1:6),
                   paste0("G", 8:24),
                   paste0("I", 1:6),
                   paste0("I", 8:24),
                   paste0("L", c(1,3:24)),
                   paste0("N", 1:24))
names(hek293t_cells) <- rep("HEK293T", length(hek293t_cells))

# filter to hek cells
storm_kb_res_gene_counts.remap.hek <- lapply(storm_kb_res_gene_counts.remap, function(x) {
  return(x[rownames(x) %in% hek293t_cells,])
})

saveRDS(storm_kb_res_gene_counts.remap.hek,
        file = "~/Documents/manuscripts/storm_seq/subsample_analysis/complexity_curves/storm/storm_hek_complexity_curve_to_plot.rds")

## plot preliminary results
library(ggplot2)
storm_kb_res_gene_counts.remap.hek <- lapply(names(storm_kb_res_gene_counts.remap.hek),
                                             function(x) {
                                               x.counts <- storm_kb_res_gene_counts.remap.hek[[x]]
                                               x.counts$depth <- x
                                               return(x.counts)
                                             })

storm_kb_res_gene_counts.remap.hek.long <- do.call(rbind,
                                                   storm_kb_res_gene_counts.remap.hek)

# filter out low expressing cells
storm_kb_res_gene_counts.remap.hek.long <- storm_kb_res_gene_counts.remap.hek.long[storm_kb_res_gene_counts.remap.hek.long$detected_genes > 200,]
storm_kb_res_gene_counts.remap.hek.long$depth <- factor(storm_kb_res_gene_counts.remap.hek.long$depth,
                                                        levels = c("50k",
                                                                   "100k",
                                                                   "150k",
                                                                   "200k",
                                                                   "250k",
                                                                   "500k",
                                                                   "1M"))

storm_kb_res_gene_counts.remap.hek.long$technology <- "STORM-seq"
ggplot(storm_kb_res_gene_counts.remap.hek.long,
       aes(y = detected_txps, x = detected_genes)) +
  geom_point(aes(color = depth), size = 3, alpha = 0.75) +
  ylab("Detected Transcripts") +
  xlab("Detected Genes") +
  scale_color_viridis_d() +
  ggtitle("STORM-seq Gene and Transcript Detection") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

## looks good, let's move on

# extract for vasa
vasa_kb_res <- gather_kb_res(getwd())
names(vasa_kb_res) <- list.dirs(full.names = FALSE,
                                 recursive = FALSE)
saveRDS(vasa_kb_res, file = "vasa_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds") 
vasa_kb_res <- readRDS("vasa/vasa_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")

# now need to go through and pull out the gene and txp counts that are not erccs
vasa_kb_res_gene_counts <- lapply(vasa_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  txp_counts <- x[["txp"]]
  txp_counts <- counts(txp_counts[grep("ERCC-",
                                       invert = TRUE,
                                       rownames(txp_counts)),])
  detected_genes <- colSums2(gene_counts >= 1)
  umi_frags <- colSums2(gene_counts)
  detected_txps <- colSums2(txp_counts >= 1)
  # assembl
  return(data.frame(detected_genes = detected_genes,
                    detected_txps = detected_txps,
                    umi_frags = umi_frags,
                    row.names = cells))
})
saveRDS(vasa_kb_res_gene_counts, file = "vasa_gene_txp_counts_ens101_with_erccs_for_plotting.rds") 
vasa_kb_res_gene_counts <- readRDS("vasa/vasa_gene_txp_counts_ens101_with_erccs_for_plotting.rds")

p1_neg_controls <- c("vai_vasa_plate1_merged_009", "vai_vasa_plate1_merged_047",
                     "vai_vasa_plate1_merged_108", "vai_vasa_plate1_merged_123",
                     "vai_vasa_plate1_merged_207", "vai_vasa_plate1_merged_318",
                     "vai_vasa_plate1_merged_334", "vai_vasa_plate1_merged_379")
p2_neg_controls <- c("vai_vasa_plate2_merged_009", "vai_vasa_plate2_merged_047",
                     "vai_vasa_plate2_merged_108", "vai_vasa_plate2_merged_123",
                     "vai_vasa_plate2_merged_207", "vai_vasa_plate2_merged_318",
                     "vai_vasa_plate2_merged_334", "vai_vasa_plate2_merged_379")

# filter the negative controls and to HEK293T cells
format_nums <- sprintf("%03d", 0:384)
vasa_plate1 <- paste0("vai_vasa_plate1_merged_",
                      format_nums)
names(vasa_plate1) <- c(rep("K562", 192),
                        rep("HEK293T", 192))
vasa_plate1 <- vasa_plate1[!vasa_plate1 %in% p1_neg_controls]
vasa_plate1.hek <- vasa_plate1[names(vasa_plate1) %in% "HEK293T"]

vasa_plate2 <- paste0("vai_vasa_plate2_merged_",
                      format_nums)
names(vasa_plate2) <- c(rep("K562", 192),
                        rep("HEK293T", 192))
vasa_plate2 <- vasa_plate2[!vasa_plate2 %in% p2_neg_controls]
vasa_plate2.hek <- vasa_plate2[names(vasa_plate2) %in% "HEK293T"]

vasa_kb_res_gene_counts.hek <- lapply(vasa_kb_res_gene_counts, function(x) {
  return(x[rownames(x) %in% c(vasa_plate1.hek,
                              vasa_plate2.hek),])
})


vasa_kb_res_gene_counts.hek <- lapply(names(vasa_kb_res_gene_counts.hek),
                                      function(x) {
                                               x.counts <- vasa_kb_res_gene_counts.hek[[x]]
                                               x.counts$depth <- x
                                               x.counts$technology <- "VASA-seq"
                                               return(x.counts)
                                             })
vasa_kb_res_gene_counts.hek.long <- do.call(rbind,
                                            vasa_kb_res_gene_counts.hek)

ggplot(vasa_kb_res_gene_counts.hek.long,
       aes(y = detected_txps, x = detected_genes)) +
  geom_point(aes(color = depth)) +
  geom_smooth(method = "loess",
              se = FALSE,
              linewidth = 2,
              color = "black") +
  scale_color_viridis_d() +
  theme_bw(12)

# extract for sstotal
sstotal_kb_res <- gather_kb_res(getwd())
names(sstotal_kb_res) <- list.dirs(full.names = FALSE,
                                   recursive = FALSE)
saveRDS(sstotal_kb_res, file = "sstotal_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds") 
sstotal_kb_res <- readRDS("sstotal/sstotal_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")

# now need to go through and pull out the gene and txp counts that are not erccs
sstotal_kb_res_gene_counts <- lapply(sstotal_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  txp_counts <- x[["txp"]]
  txp_counts <- counts(txp_counts[grep("ERCC-",
                                       invert = TRUE,
                                       rownames(txp_counts)),])
  detected_genes <- colSums2(gene_counts >= 1)
  umi_frags <- colSums2(gene_counts)
  detected_txps <- colSums2(txp_counts >= 1)
  # assemble
  return(data.frame(detected_genes = detected_genes,
                    detected_txps = detected_txps,
                    umi_frags = umi_frags,
                    row.names = cells))
})
saveRDS(sstotal_kb_res_gene_counts, file = "sstotal_gene_txp_counts_ens101_with_erccs_for_plotting.rds") 
sstotal_kb_res_gene_counts <- readRDS("sstotal/sstotal_gene_txp_counts_ens101_with_erccs_for_plotting.rds")
sstotal_kb_res_gene_counts.hek <- lapply(names(sstotal_kb_res_gene_counts),
                                         function(x) {
                                           x.counts <- sstotal_kb_res_gene_counts[[x]]
                                           x.counts$depth <- x
                                           x.counts$technology <- "Smart-seq-total"
                                           return(x.counts)
                                         })
sstotal_kb_res_gene_counts.hek.long <- do.call(rbind,
                                               sstotal_kb_res_gene_counts.hek)

# filter out low expressing cells
sstotal_kb_res_gene_counts.hek.long <- sstotal_kb_res_gene_counts.hek.long[sstotal_kb_res_gene_counts.hek.long$detected_genes > 100,]

# extract for ss3xpress
ss3xpress_kb_res <- gather_kb_res(getwd())
names(ss3xpress_kb_res) <- list.dirs(full.names = FALSE,
                                   recursive = FALSE)
saveRDS(ss3xpress_kb_res, file = "ss3xpress_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds") 
ss3xpress_kb_res <- readRDS("ss3x/ss3xpress_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")

# now need to go through and pull out the gene and txp counts that are not erccs
# there are not ERCCs in these data, but a carryover from above
ss3xpress_kb_res_gene_counts <- lapply(ss3xpress_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  txp_counts <- x[["txp"]]
  txp_counts <- counts(txp_counts[grep("ERCC-",
                                       invert = TRUE,
                                       rownames(txp_counts)),])
  detected_genes <- colSums2(gene_counts >= 1)
  umi_frags <- colSums2(gene_counts)
  detected_txps <- colSums2(txp_counts >= 1)
  # assemble
  return(data.frame(detected_genes = detected_genes,
                    detected_txps = detected_txps,
                    umi_frags = umi_frags,
                    row.names = cells))
})
saveRDS(ss3xpress_kb_res_gene_counts, file = "ss3xpress_gene_txp_counts_ens101_with_erccs_for_plotting.rds") 
ss3xpress_kb_res_gene_counts <- readRDS("ss3x/ss3xpress_gene_txp_counts_ens101_with_erccs_for_plotting.rds")
ss3xpress_kb_res_gene_counts.hek <- lapply(names(ss3xpress_kb_res_gene_counts),
                                         function(x) {
                                           x.counts <- ss3xpress_kb_res_gene_counts[[x]]
                                           x.counts$depth <- x
                                           x.counts$technology <- "Smart-seq3xpress"
                                           return(x.counts)
                                         })
ss3xpress_kb_res_gene_counts.hek.long <- do.call(rbind,
                                               ss3xpress_kb_res_gene_counts.hek)

# combine and plot
all_tech <- rbind(storm_kb_res_gene_counts.remap.hek.long,
                  vasa_kb_res_gene_counts.hek.long,
                  sstotal_kb_res_gene_counts.hek.long,
                  ss3xpress_kb_res_gene_counts.hek.long)
all_tech$depth <- factor(all_tech$depth,
                         levels = c("50k",
                                    "100k",
                                    "150k",
                                    "200k",
                                    "250k",
                                    "500k",
                                    "1M"))

# truncate at 150k for comparisons due to sampling issues
all_tech <- all_tech[all_tech$depth %in% c("50k",
                                           "100k",
                                           "150k"),]

# re-level for plotting
ggplot(all_tech,
       aes(y = detected_genes, x = depth,
           color = technology, group = technology)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymax = ..y.., ymin = ..y.., group = technology),
    width = 0.75,
    position = position_dodge(width = 0.75),
    size = 2,
    color = "black"
  ) +
  ylab("Detected Genes") +
  xlab("Read Depth") +
  ggtitle("Gene Detection") +
  scale_y_continuous(limits = c(0, 12000), breaks = seq(0,12000,3000)) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                               `VASA-seq`="#CF5917FF",
                               `Smart-seq-total` = "#8FC2FD",
                               `Smart-seq3xpress` = "#8CCC98"),
                    name = "Technology") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

ggplot(all_tech,
       aes(y = detected_txps, x = depth,
           color = technology, group = technology)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymax = ..y.., ymin = ..y.., group = technology),
    width = 0.75,
    position = position_dodge(width = 0.75),
    size = 2,
    color = "black"
  ) +
  ylab("Detected Transcripts") +
  xlab("Read Depth") +
  ggtitle("Transcript Detection") +
  scale_y_continuous(limits = c(0, 16000), breaks = seq(0,16000,4000)) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF",
                                `Smart-seq-total` = "#8FC2FD",
                                `Smart-seq3xpress` = "#8CCC98"),
                                name = "Technology") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5)
  )