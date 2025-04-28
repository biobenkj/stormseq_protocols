## biotype diversity plots
## import the GTF
library(rtracklayer)

gtf <- rtracklayer::import("Homo_sapiens.GRCh38.101.gtf")
gtf.gene <- gtf[gtf$type %in% "gene",]
gtf.gene.df <- data.frame(gene_id = gtf.gene$gene_id,
                          gene_biotype = gtf.gene$gene_biotype)

gene_biotypes <- gtf.gene.df %>% select(gene_id, gene_biotype)
biotype_mapping <- function(biotype) {
  if (biotype %in% c("protein_coding")) {
    return("protein_coding")
  } else if (biotype %in% c("lncRNA")) {
    return("lncRNA")
  } else if (biotype %in% c("miscRNA")) {
    return("miscRNA")
  } else if (biotype %in% c("pseudogene")) {
    return("pseudogenes")
  } else if (biotype %in% c("rRNA", "MT_rRNA")) {
    return("rRNA")
  } else if (biotype %in% c("miRNA")) {
    return("miRNA")
  } else if (biotype %in% c("snoRNA")) {
    return("snoRNA")
  } else if (biotype %in% c("sncRNA")) {
    return("sncRNA")
  } else if (biotype %in% c("snRNA")) {
    return("snRNA")
  } else if (biotype %in% c("scaRNA")) {
    return("scaRNA")
  } else if (biotype %in% c("sRNA")) {
    return("sRNA")
  } else {
    return("other")
  }
}
gtf.gene.df$mapped_biotype <- sapply(gtf.gene.df$gene_biotype, biotype_mapping)

# per reviewer request, restrict to protein coding
gtf.gene.df.pc <- gtf.gene.df[gtf.gene.df$mapped_biotype %in% "protein_coding",]

# loop back through each technology and restrict to protein coding genes
# STORM
storm_kb_res <- readRDS("storm/storm_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
storm_kb_res_gene_counts.pc <- lapply(storm_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  gene_counts <- fix_rownames(gene_counts)
  # restrict to protein_coding
  gene_counts <- gene_counts[rownames(gene_counts) %in% gtf.gene.df.pc$gene_id,]
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

storm_kb_res_gene_counts.pc.remap <- lapply(storm_kb_res_gene_counts.pc, function(x) {
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
storm_kb_res_gene_counts.pc.remap.hek <- lapply(storm_kb_res_gene_counts.pc.remap, function(x) {
  return(x[rownames(x) %in% hek293t_cells,])
})

storm_kb_res_gene_counts.pc.remap.hek <- lapply(names(storm_kb_res_gene_counts.pc.remap.hek),
                                                function(x) {
                                                  x.counts <- storm_kb_res_gene_counts.pc.remap.hek[[x]]
                                                  x.counts$depth <- x
                                                  return(x.counts)
                                                })

storm_kb_res_gene_counts.pc.remap.hek.long <- do.call(rbind,
                                                      storm_kb_res_gene_counts.pc.remap.hek)

# filter out low expressing cells
storm_kb_res_gene_counts.pc.remap.hek.long <- storm_kb_res_gene_counts.pc.remap.hek.long[storm_kb_res_gene_counts.pc.remap.hek.long$detected_genes > 200,]
storm_kb_res_gene_counts.pc.remap.hek.long$depth <- factor(storm_kb_res_gene_counts.pc.remap.hek.long$depth,
                                                           levels = c("50k",
                                                                      "100k",
                                                                      "150k",
                                                                      "200k",
                                                                      "250k",
                                                                      "500k",
                                                                      "1M"))

storm_kb_res_gene_counts.pc.remap.hek.long$technology <- "STORM-seq"
ggplot(storm_kb_res_gene_counts.pc.remap.hek.long,
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

# per reviewer request, look at all gene biotypes found across the different cell types
# first, sum up the gene counts across cells for detection rate
# remap the wells to the proper place and remove controls
storm_remap <- read.delim("storm/well_map.txt",
                          header = FALSE)
storm_kb_res.remap <- lapply(storm_kb_res, function(x) {
  # grab the genes
  x.gene <- x$gene
  colnames(x.gene) <- gsub("\\/", "", colnames(x.gene))
  x.gene <- fix_rownames(x.gene)
  storm_remap.match <- match(colnames(x.gene),
                             storm_remap$V1)
  storm_remap.tmp <- storm_remap[storm_remap.match,]
  stopifnot(all(storm_remap.tmp$V1 == colnames(x.gene)))
  colnames(x.gene) <- storm_remap.tmp$V2
  return(x.gene)
})


# annotate cell types
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

# filter out controls
storm_kb_res.remap <- lapply(storm_kb_res.remap, function(x) {
  return(x[,!colnames(x) %in% c(pos_control,
                                neg_control)])
})

all_celltypes <- c(k562_cells,
                   hek293t_cells,
                   rmg2_cells)

# organize
storm_kb_res.remap <- lapply(storm_kb_res.remap, function(x) {
  x.match <- match(colnames(x),
                   all_celltypes)
  all_celltypes.tmp <- all_celltypes[x.match]
  stopifnot(all(all_celltypes.tmp == colnames(x)))
  x$cell_type <- names(all_celltypes.tmp)
  return(x)
})

# flatten
storm_kb_res.remap.flat <- lapply(storm_kb_res.remap, function(x) {
  # split by cell type
  hek <- x[,x$cell_type %in% "HEK293T"]
  k562 <- x[,x$cell_type %in% "K562"]
  rmg2 <- x[,x$cell_type %in% "RMG2"]
  # calculate detection rate
  hek.detect <- hek[rowSums2(counts(hek)) > 0,]
  k562.detect <- k562[rowSums2(counts(k562)) > 0,]
  rmg2.detect <- rmg2[rowSums2(counts(rmg2)) > 0,]
  # flatten
  rowData(hek.detect)$cell_type <- "HEK293T"
  rowData(k562.detect)$cell_type <- "K562"
  rowData(rmg2.detect)$cell_type <- "RMG2"
  detect.df <- data.frame(gene_ids = c(rownames(hek.detect),
                                       rownames(k562.detect),
                                       rownames(rmg2.detect)),
                          cell_type = c(rowData(hek.detect)$cell_type,
                                        rowData(k562.detect)$cell_type,
                                        rowData(rmg2.detect)$cell_type))
  return(detect.df)
})

storm_kb_res.remap.flat.biotypes <- lapply(storm_kb_res.remap.flat, function(x) {
  # restrict to what is found in the GTF biotypes
  x <- x[x$gene_ids %in% gtf.gene.df$gene_id,]
  # restrict to all biotypes mapped
  biotypes_detected <- do.call(c, lapply(x$gene_ids, function(y) {
    gtf.gene.df.sub <- gtf.gene.df[gtf.gene.df$gene_id %in% y,]
    return(gtf.gene.df.sub$mapped_biotype)
  }))
  # assemble
  x$mapped_biotype <- biotypes_detected
  message("Done.")
  return(x)
})

storm_kb_res.remap.flat.biotypes.names <- lapply(names(storm_kb_res.remap.flat.biotypes),
                                                 function(x) {
                                                   x.counts <- storm_kb_res.remap.flat.biotypes[[x]]
                                                   x.counts$depth <- x
                                                   return(x.counts)
                                                 })

storm_kb_res.remap.flat.biotypes.names.long <- do.call(rbind,
                                                       storm_kb_res.remap.flat.biotypes.names)

# factor depth for faceting
storm_kb_res.remap.flat.biotypes.names.long$depth <- factor(storm_kb_res.remap.flat.biotypes.names.long$depth,
                                                            levels = c("50k",
                                                                       "100k",
                                                                       "150k",
                                                                       "200k",
                                                                       "250k",
                                                                       "500k",
                                                                       "1M"))

# summarize
# Summarize counts
df_summary <- storm_kb_res.remap.flat.biotypes.names.long %>%
  group_by(cell_type, depth, mapped_biotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_type, depth) %>%
  mutate(percent = (count / sum(count)) * 100)  # Convert to percentage

df_summary$mapped_biotype <- ifelse(df_summary$mapped_biotype %in% "other",
                                    "multi-annotated",
                                    df_summary$mapped_biotype) 

# order the factors
df_summary$mapped_biotype <- factor(df_summary$mapped_biotype,
                                    levels = c("protein_coding",
                                               "lncRNA",
                                               "multi-annotated",
                                               "miRNA",
                                               "pseudogenes",
                                               "rRNA",
                                               "scaRNA",
                                               "snoRNA",
                                               "snRNA"))

# Plot
ggplot(df_summary, aes(x = depth, y = percent, color = mapped_biotype, group = mapped_biotype)) +
  geom_point(size = 3) +  # Dot plot
  geom_line(linewidth = 1) +   # Connecting lines
  facet_wrap(~cell_type) +
  scale_color_viridis_d(option = "plasma", 
                        end = 0.85,
                        labels = c("protein_coding" = "Protein coding",
                                   "lncRNA" = "lncRNA",
                                   "multi-annotated" = "Multi-annotated",
                                   "miRNA" = "miRNA",
                                   "pseudogenes" = "Pseudogenes",
                                   "rRNA" = "rRNA",
                                   "scaRNA" = "scaRNA",
                                   "snoRNA" = "snoRNA",
                                   "snRNA" = "snRNA")) +
  labs(
    title = "Mapped Biotype Proportion Across Depths",
    x = "Depth",
    y = "Percentage",
    color = "Mapped Biotype"
  ) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(storm_kb_res_gene_counts.pc.remap.hek.long,
       aes(y = detected_txps, x = detected_genes)) +
  geom_point(aes(color = depth), size = 3, alpha = 0.75) +
  # geom_smooth(method = "loess",
  #             se = FALSE,
  #             linewidth = 2,
  #             color = "black") +
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

# VASA
vasa_kb_res <- readRDS("vasa/vasa_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
vasa_kb_res_gene_counts.pc <- lapply(vasa_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  gene_counts <- fix_rownames(gene_counts)
  # restrict to protein_coding
  gene_counts <- gene_counts[rownames(gene_counts) %in% gtf.gene.df.pc$gene_id,]
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

vasa_kb_res_gene_counts.pc.hek <- lapply(vasa_kb_res_gene_counts.pc, function(x) {
  return(x[rownames(x) %in% c(vasa_plate1.hek,
                              vasa_plate2.hek),])
})

vasa_kb_res_gene_counts.pc.hek <- lapply(names(vasa_kb_res_gene_counts.pc.hek),
                                         function(x) {
                                           x.counts <- vasa_kb_res_gene_counts.pc.hek[[x]]
                                           x.counts$depth <- x
                                           x.counts$technology <- "VASA-seq"
                                           return(x.counts)
                                         })
vasa_kb_res_gene_counts.pc.hek.long <- do.call(rbind,
                                               vasa_kb_res_gene_counts.pc.hek)

ggplot(vasa_kb_res_gene_counts.pc.hek.long,
       aes(y = detected_txps, x = detected_genes)) +
  geom_point(aes(color = depth)) +
  geom_smooth(method = "loess",
              se = FALSE,
              linewidth = 2,
              color = "black") +
  scale_color_viridis_d() +
  theme_bw(12)

# SStotal
sstotal_kb_res <- readRDS("sstotal/sstotal_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
sstotal_kb_res_gene_counts.pc <- lapply(sstotal_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  gene_counts <- fix_rownames(gene_counts)
  # restrict to protein_coding
  gene_counts <- gene_counts[rownames(gene_counts) %in% gtf.gene.df.pc$gene_id,]
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

sstotal_kb_res_gene_counts.pc.hek <- lapply(names(sstotal_kb_res_gene_counts.pc),
                                            function(x) {
                                              x.counts <- sstotal_kb_res_gene_counts.pc[[x]]
                                              x.counts$depth <- x
                                              x.counts$technology <- "Smart-seq-total"
                                              return(x.counts)
                                            })
sstotal_kb_res_gene_counts.pc.hek.long <- do.call(rbind,
                                                  sstotal_kb_res_gene_counts.pc.hek)

# filter out low expressing cells
sstotal_kb_res_gene_counts.pc.hek.long <- sstotal_kb_res_gene_counts.pc.hek.long[sstotal_kb_res_gene_counts.pc.hek.long$detected_genes > 100,]

# SS3x
ss3xpress_kb_res <- readRDS("ss3x/ss3xpress_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
ss3xpress_kb_res_gene_counts.pc <- lapply(ss3xpress_kb_res, function(x) {
  gene_counts <- x[["gene"]]
  cells <- colnames(gene_counts)
  cells <- gsub("\\/", "", cells)
  gene_counts <- counts(gene_counts[grep("ERCC-",
                                         invert = TRUE,
                                         rownames(gene_counts)),])
  gene_counts <- fix_rownames(gene_counts)
  # restrict to protein_coding
  gene_counts <- gene_counts[rownames(gene_counts) %in% gtf.gene.df.pc$gene_id,]
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

ss3xpress_kb_res_gene_counts.pc.hek <- lapply(names(ss3xpress_kb_res_gene_counts.pc),
                                              function(x) {
                                                x.counts <- ss3xpress_kb_res_gene_counts.pc[[x]]
                                                x.counts$depth <- x
                                                x.counts$technology <- "Smart-seq3xpress"
                                                return(x.counts)
                                              })
ss3xpress_kb_res_gene_counts.pc.hek.long <- do.call(rbind,
                                                    ss3xpress_kb_res_gene_counts.pc.hek)

# combine and plot
all_tech <- rbind(storm_kb_res_gene_counts.pc.remap.hek.long,
                  vasa_kb_res_gene_counts.pc.hek.long,
                  sstotal_kb_res_gene_counts.pc.hek.long,
                  ss3xpress_kb_res_gene_counts.pc.hek.long)
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
  ylab("Detected Protein-coding Genes") +
  xlab("Read Depth") +
  ggtitle("Protein-coding Gene Detection") +
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