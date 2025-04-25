## process the TEProf3 output

# pull in the bulk total RNA short-read data that has been run in guided mode
# using the STORM K-562 pseudobulk consensus TE GTF
bulk_tes_storm_consensus <- read.delim("bulk/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")

bulk_tes_storm_consensus_quant <- read.delim("bulk/teprof3_output_quantification.TE.tsv.gz")

# filter based on recommended TEProf3 docs
library(dplyr)

bulk_tes_storm_consensus_quant.filt <- bulk_tes_storm_consensus_quant %>%
  group_by(transcript_id) %>%
  filter(
    all(stringtie_tpm >= 1),
    all(perfect_SJ_uniqlymapped_read_downstream >= 1),
    all(SJ_uniqlymapped_read_upstream <= (0.5 * perfect_SJ_uniqlymapped_read_downstream))
  ) %>%
  ungroup()

table(length(unique(bulk_tes_storm_consensus_quant.filt$transcript_id)))
# 440

# quick look
plot(hist(bulk_tes_storm_consensus_quant.filt$sum_tpm))

## STORM
storm_k562_sc_tes <- read.delim("storm/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")
storm_k562_sc_tes.quant <- read.delim("storm/teprof3_output_quantification.TE.tsv.gz")

# pull in the expression SCE from the benchmarking to filter to QC'd cells
storm_kb_res <- readRDS("storm/storm_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
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

# use 1M subsample
storm_kb_res.remap <- storm_kb_res.remap$`1M`

# pull out the ERCCs
ercc.sce <- storm_kb_res.remap[grep("ERCC-", rownames(storm_kb_res.remap)),]
altExp(storm_kb_res.remap, "ERCC") <- ercc.sce
storm_kb_res.remap <- storm_kb_res.remap[grep("ERCC-", rownames(storm_kb_res.remap), invert = TRUE),]

# re-calculate num genes
storm_kb_res.remap$NumGenesExpressed <- colSums2(counts(storm_kb_res.remap) > 0)

# QC
library(scuttle)
storm_kb_res.remap <- addPerCellQCMetrics(storm_kb_res.remap)
filt_reasons <- perCellQCFilters(storm_kb_res.remap,
                                 sub.fields = c("altexps_ERCC_percent"))
colSums(as.matrix(filt_reasons))
# low_lib_size            low_n_features high_altexps_ERCC_percent                   discard 
# 2                         4                         8                        10

# filter
storm_sce.filt <- storm_kb_res.remap[,!filt_reasons$discard]
table(storm_sce.filt$cell_type)
# HEK293T    K562    RMG2 
# 128      96      97 

storm_sce.filt <- storm_sce.filt[,storm_sce.filt$cell_type %in% "K562"]

qcd_k562 <- paste0(colnames(storm_sce.filt),
                   "_teprof3")

# need to remap the TEProf3 storm te column names....
storm_remap <- read.delim("storm/well_map.txt",
                          header = FALSE)
storm_remap$V1 <- paste0(storm_remap$V1,
                         "_teprof3")
storm_remap$V2 <- paste0(storm_remap$V2,
                         "_teprof3")

# remap in a vectorized fashion
storm_remap <- setNames(storm_remap$V2, storm_remap$V1)

storm_k562_sc_tes.quant$remapped <- storm_remap[storm_k562_sc_tes.quant$sample]

# filter
storm_k562_sc_tes.quant <- storm_k562_sc_tes.quant[storm_k562_sc_tes.quant$remapped %in% qcd_k562,]

# process and filter
library(dplyr)

# Count the number of non-zero `sum_tpm` values per transcript_id
# be inclusive for comparison to bulk
storm_k562_sc_tes.quant.filt <- storm_k562_sc_tes.quant %>%
  group_by(transcript_id) %>%
  mutate(total_cells = n(),  # Count total occurrences of each transcript_id
         nonzero_cells = sum(stringtie_tpm > 0)) %>%  # Count non-zero `sum_tpm` values
  mutate(perfect_sj_nonzero_cells = sum(perfect_SJ_uniqlymapped_read_downstream > 0)) %>%
  filter(nonzero_cells / total_cells >= 0.1) %>%  # Keep if at least 10% have non-zero sum_tpm
  ungroup()
length(unique(storm_k562_sc_tes.quant.filt$transcript_id))
# 663


# let's just subset to bulk TE-txps since that's our "ground truth"
storm_single_cell_oncoexaptations.bulk <- storm_k562_sc_tes.quant[storm_k562_sc_tes.quant$transcript_id %in% bulk_tes_storm_consensus_quant.filt$transcript_id,]


# alright... let's tease these apart
library(tidyr)
# aggregate stringtie_tpm per gene_name and sample
# I'm from the future (below code) to re-filter to desired TE-txps prior to collapsing to gene!
# storm_single_cell_oncoexaptations.bulk.filt <- storm_single_cell_oncoexaptations.bulk[storm_single_cell_oncoexaptations.bulk$transcript_id %in% storm_k562_sc_tes.common_bulk.tes$transcript_id,]
# storm_filt_txp_ids <- storm_k562_sc_tes.common_bulk.tes$transcript_id
# names(storm_filt_txp_ids) <- storm_k562_sc_tes.common_bulk.tes$te_name_to_plot

# filter out mono exons too!
mono_exons <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/storm_pseudobulk/mono_exons_with_header.txt")
# table(mono_exons$transcript_id %in% storm_single_cell_oncoexaptations.bulk.filt$transcript_id)
# FALSE 
# 12930
# good.

# discovery and curated
storm_k562_sc_tes_aggregated_data <- storm_single_cell_oncoexaptations.bulk %>%
  group_by(gene_name, remapped) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")  # Sum over all transcript_id

# DIFFERENT AFTER DEDUP!!!!
# storm_k562_sc_tes_aggregated_data <- storm_single_cell_oncoexaptations.bulk.filt %>%
#   group_by(transcript_id, remapped) %>%
#   summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")

# make the expression matrix
storm_k562_sc_tes_expression_matrix <- storm_k562_sc_tes_aggregated_data %>%
  spread(key=remapped, value=sum_tpm, fill=0)  # Convert to wide format

# AFTER DEDUP!
# ens_id.m <- match(storm_k562_sc_tes_expression_matrix$transcript_id,
#                   storm_filt_txp_ids)
# storm_filt_txp_ids <- storm_filt_txp_ids[ens_id.m]
# storm_k562_sc_tes_expression_matrix$transcript_id <- names(storm_filt_txp_ids)
# 
# storm_k562_sc_tes_expression_matrix <- cbind(storm_k562_sc_tes_expression_matrix$transcript_id,
#                                              log2(storm_k562_sc_tes_expression_matrix[,c(2:97)] + 1))

storm_k562_sc_tes_expression_matrix <- cbind(storm_k562_sc_tes_expression_matrix$gene_name,
                                             log2(storm_k562_sc_tes_expression_matrix[,c(2:97)] + 1))
colnames(storm_k562_sc_tes_expression_matrix)[1] <- "gene_name"

# Exclude any cell with less than 9k genes
storm_gene_counts <- storm_sce.filt$NumGenesExpressed
names(storm_gene_counts) <- paste0(colnames(storm_sce.filt),
                                   "_teprof3")
storm_gene_counts <- storm_gene_counts[storm_gene_counts > 9000]
storm_k562_sc_tes_expression_matrix <- storm_k562_sc_tes_expression_matrix[,colnames(storm_k562_sc_tes_expression_matrix) %in% c(names(storm_gene_counts),
                                                                                                                                 "gene_name")]

# storm_k562_sc_tes_all_expression_matrix.log <- storm_k562_sc_tes.quant.filt_all %>%
#   select(transcript_id, remapped, stringtie_tpm) %>%
#   mutate(log2_tpm = log2(stringtie_tpm + 1)) %>%
#   pivot_wider(id_cols = transcript_id,  # Use only transcript_id as identifier
#               names_from = remapped,
#               values_from = log2_tpm,
#               values_fill = list(log2_tpm = 0))
# storm_k562_sc_tes_all_expression_matrix.log <- as.data.frame(storm_k562_sc_tes_all_expression_matrix.log)
# rownames(storm_k562_sc_tes_all_expression_matrix.log) <- storm_k562_sc_tes_all_expression_matrix.log$transcript_id
# storm_k562_sc_tes_all_expression_matrix.log$transcript_id <- NULL

# do the same for the bulk
# pull out all the common TE-derived transcripts from bulk
# bulk_tes_storm_consensus_quant.filt.common_tu <- bulk_tes_storm_consensus_quant.filt[bulk_tes_storm_consensus_quant.filt$transcript_id %in% storm_all_cell_oncoexaptations,]

# aggregate stringtie_tpm per gene_name and sample
# I'm from the future (below code) to re-filter to desired TE-txps prior to collapsing to gene!
# bulk_tes_storm_consensus_quant.filt.tefilt <- bulk_tes_storm_consensus_quant.filt[bulk_tes_storm_consensus_quant.filt$transcript_id %in% storm_k562_sc_tes.common_bulk.tes$transcript_id,]
# bulk_filt_txp_ids <- storm_k562_sc_tes.common_bulk.tes$transcript_id
# names(bulk_filt_txp_ids) <- storm_k562_sc_tes.common_bulk.tes$te_name_to_plot

bulk_k562_tes_aggregated_data <- bulk_tes_storm_consensus_quant.filt %>%
  group_by(gene_name, sample) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")  # Sum over all transcript_id

## DIFFERENT AFTER DEDUP!!!!
# bulk_k562_tes_aggregated_data <- bulk_tes_storm_consensus_quant.filt.tefilt %>%
#   group_by(transcript_id, sample) %>%
#   summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")

# make the expression matrix
bulk_k562_tes_expression_matrix <- bulk_k562_tes_aggregated_data %>%
  spread(key=sample, value=sum_tpm, fill=0)  # Convert to wide format

# AFTER DEDUP!
# ens_id.m <- match(bulk_k562_tes_expression_matrix$transcript_id,
#                   bulk_filt_txp_ids)
# bulk_filt_txp_ids <- bulk_filt_txp_ids[ens_id.m]
# bulk_k562_tes_expression_matrix$transcript_id <- names(bulk_filt_txp_ids)
# bulk_k562_tes_expression_matrix <- cbind(bulk_k562_tes_expression_matrix$transcript_id,
#                                          log2(bulk_k562_tes_expression_matrix[,c(2:3)] + 1))

bulk_k562_tes_expression_matrix <- cbind(bulk_k562_tes_expression_matrix$gene_name,
                                         log2(bulk_k562_tes_expression_matrix[,c(2:3)] + 1))
colnames(bulk_k562_tes_expression_matrix)[1] <- "gene_name"

# bulk_k562_sc_tes_all_expression_matrix.log <- bulk_tes_storm_consensus_quant.filt.common_tu %>%
#   select(transcript_id, sample, stringtie_tpm) %>%
#   mutate(log2_tpm = log2(stringtie_tpm + 1)) %>%
#   pivot_wider(id_cols = transcript_id,  # Use only transcript_id as identifier
#               names_from = sample,
#               values_from = log2_tpm,
#               values_fill = list(log2_tpm = 0))
# bulk_k562_sc_tes_all_expression_matrix.log <- as.data.frame(bulk_k562_sc_tes_all_expression_matrix.log)
# rownames(bulk_k562_sc_tes_all_expression_matrix.log) <- bulk_k562_sc_tes_all_expression_matrix.log$transcript_id
# bulk_k562_sc_tes_all_expression_matrix.log$transcript_id <- NULL

#bulk_and_storm_sc_tes <- inner_join(bulk_k562_tes_expression_matrix, storm_k562_sc_tes_expression_matrix, by = "gene_name")
bulk_and_storm_sc_tes <- inner_join(bulk_k562_tes_expression_matrix,
                                    storm_k562_sc_tes_expression_matrix, by = "gene_name")
rownames(bulk_and_storm_sc_tes) <- bulk_and_storm_sc_tes$gene_name
bulk_and_storm_sc_tes$gene_name <- NULL
# bulk_and_storm_sc_tes$storm_pseudobulk <- rowMeans(bulk_and_storm_sc_tes[,c(4:99)])

# go through and find the TE-Gene combo for annotations
# e.g. L2c-ENSG00000229140
# pull in the GTF and find out what genes we have going on
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.101.ercc92patched.gtf.gz")
gtf.genes <- gtf[mcols(gtf)$type == "gene",]

storm_k562_sc_tes.common_bulk <- storm_k562_sc_tes[storm_k562_sc_tes$gene_name %in% rownames(bulk_and_storm_sc_tes),]
storm_k562_sc_tes.common_bulk <- storm_k562_sc_tes.common_bulk[storm_k562_sc_tes.common_bulk$transcript_id %in% storm_single_cell_oncoexaptations.bulk$transcript_id,]
# deduplicate - see notes below on duplicate TE transcripts
# also reinforces the use of `sum_tpm` instead of transcript level
# this should probably be addressed in TEProf3 proper at some point
storm_k562_sc_tes.common_bulk_dedupe <- storm_k562_sc_tes.common_bulk %>%
  group_by(gene_name) %>%
  group_modify(~ {
    if (n_distinct(.x$TE_subfamily) == 1) {
      # If TE_subfamily is identical, collapse to one row (first row)
      .x[1, ]
    } else {
      # Otherwise, return all rows for this tss_gene_name
      .x
    }
  }) %>%
  ungroup()

# finer grained dupe resolution
storm_k562_sc_tes.common_bulk_dedupe.fine <- storm_k562_sc_tes.common_bulk_dedupe[duplicated(storm_k562_sc_tes.common_bulk_dedupe$gene_name),]
length(unique(storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name))
# 6
# ENSG00000138190
# ENSG00000164292
# ENSG00000226419
# ENSG00000232884
# ENSG00000250337
# ENSG00000267243

duped_genes <- unique(storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name)
storm_k562_sc_tes.common_bulk_dedupe.fine <- storm_single_cell_oncoexaptations.bulk[storm_single_cell_oncoexaptations.bulk$gene_name %in% duped_genes,]

# ENSG00000103657
# ens_1 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000103657",]
# ens_1 <- ens_1[ens_1$percentage_of_expression > 0,]
# table(ens_1$transcript_id)
# # TU4690 TU4695
# # 96      9
# # TU4690 is dominant
# tu_keep <- c("TU4690")

# ENSG00000138190
ens_2 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000138190",]
ens_2 <- ens_2[ens_2$percentage_of_expression > 0,]
table(ens_2$transcript_id)
# TU348 TU354 TU360 
# 2    53    55 
tu_check <- c("TU354", "TU360")
ens_2 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
ens_2 <- ens_2[ens_2$percentage_of_expression > 0,]
table(ens_2$transcript_id)
# TU354 is dominant by expression
tu_keep <- c("TU354")

# ENSG00000164292
ens_2.1 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000164292",]
ens_2.1 <- ens_2.1[ens_2.1$percentage_of_expression > 0,]
table(ens_2.1$transcript_id)
# TU1165 TU1222 TU1226 
# 13      3     27 
tu_check <- c("TU1165", "TU1226 ")
ens_2.1 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
ens_2.1 <- ens_2.1[ens_2.1$percentage_of_expression > 0,]
table(ens_2.1$transcript_id)
# TU1165 is dominant
tu_keep <- c(tu_keep, "TU1165")

# ENSG00000226419
ens_3 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000226419",]
ens_3 <- ens_3[ens_3$percentage_of_expression > 0,]
table(ens_3$transcript_id)
# TU60 TU63 TU66 
# 6   60   18
# check again
tu_check <- unique(ens_3$transcript_id)
ens_3 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
ens_3 <- ens_3[ens_3$percentage_of_expression > 50,]
table(ens_3$transcript_id)
# TU63 is dominant
tu_keep <- c(tu_keep, "TU63")

# ENSG00000232884
ens_5 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000232884",]
ens_5 <- ens_5[ens_5$percentage_of_expression > 0,]
table(ens_5$transcript_id)
# TU918 TU922 TU924 TU948 TU955 
# 33     8    49    67    31
tu_check <- unique(ens_5$transcript_id)
ens_5 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
ens_5 <- ens_5[ens_5$percentage_of_expression > 50,]
table(ens_5$transcript_id)
# These are separate TE-derived txps that have overlapping genes...
# TU948 is an L1
# TU924 is an L2
# Keep both!
tu_keep <- c(tu_keep,
             "TU948",
             "TU924")

# ENSG00000250337
ens_7 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000250337",]
ens_7 <- ens_7[ens_7$percentage_of_expression > 0,]
table(ens_7$transcript_id)
# TU1199 TU1201 
# 47     17
tu_check <- unique(ens_7$transcript_id)
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]

# TU1199 is dominant and has a CAGE and PROcap peak
tu_keep <- c(tu_keep, "TU1199")

# ENSG00000267243
ens_8 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000267243",]
ens_8 <- ens_8[ens_8$percentage_of_expression > 0,]
table(ens_8$transcript_id)
# TU732 TU735 TU738 TU742 
# 90     3    30    22 
# double check...
tu_check <- unique(ens_8$transcript_id)
ens_8 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
ens_8 <- ens_8[ens_8$percentage_of_expression > 50,]
table(ens_8$transcript_id)

# TU742 is it's own thing
# TU732 is the same as the other 4 other than TU742
# based on IGV and CAGE/Pro-cap peaks
tu_keep <- c(tu_keep, "TU742", "TU732")

# reassemble
storm_k562_sc_tes.common_bulk_deduped <- storm_k562_sc_tes.common_bulk_dedupe %>%
  group_by(gene_name) %>%
  filter(n() == 1) %>%
  ungroup()
storm_k562_sc_tes.common_bulk_duped.fine <- storm_k562_sc_tes.common_bulk_dedupe[storm_k562_sc_tes.common_bulk_dedupe$transcript_id %in% tu_keep,]
storm_k562_sc_tes.common_bulk_deduped.clean <- rbind(storm_k562_sc_tes.common_bulk_deduped,
                                                     storm_k562_sc_tes.common_bulk_duped.fine)

# find the gene names in the gtf
gtf.genes.tes <- gtf.genes[gtf.genes$gene_id %in% storm_k562_sc_tes.common_bulk_deduped.clean$gene_name,]
gtf.genes.tes.m <- match(storm_k562_sc_tes.common_bulk_deduped.clean$gene_name,
                         gtf.genes.tes$gene_id)
gtf.genes.tes <- gtf.genes.tes[gtf.genes.tes.m,]
all(gtf.genes.tes$gene_id == storm_k562_sc_tes.common_bulk_deduped.clean$gene_name)
# TRUE
te_ens_genes <- paste0(storm_k562_sc_tes.common_bulk_deduped.clean$TE_subfamily,
                       "-",
                       gtf.genes.tes$gene_name)
storm_k562_sc_tes.common_bulk_deduped.clean$te_name_to_plot <- te_ens_genes

# pull out the final TU transcripts
storm_k562_sc_tes.common_bulk.tu <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% rownames(bulk_and_storm_sc_tes),]
storm_k562_sc_tes.common_bulk.tu$te_name_to_plot <- storm_k562_sc_tes.common_bulk.tu$teprof3_gene_name

# combine 
storm_k562_sc_tes.common_bulk.tes <- rbind(storm_k562_sc_tes.common_bulk_deduped.clean,
                                           storm_k562_sc_tes.common_bulk.tu)


## GO BACK AND RE-FILTER USING TRANSCRIPT IDS AND RECOMBINE FOR EXPRESSION!!!
# storm_k562_sc_tes.common_bulk.tes$gene_name <- ifelse(storm_k562_sc_tes.common_bulk.tes$gene_name != ".",
#                                                       storm_k562_sc_tes.common_bulk.tes$gene_name,
#                                                       storm_k562_sc_tes.common_bulk.tes$transcript_id)
# order
# storm_k562_sc_tes.common_bulk.tes.m <- match(rownames(bulk_and_storm_sc_tes),
#                                              storm_k562_sc_tes.common_bulk.tes$gene_name)
# storm_k562_sc_tes.common_bulk.tes <- storm_k562_sc_tes.common_bulk.tes[storm_k562_sc_tes.common_bulk.tes.m,]
# rownames(bulk_and_storm_sc_tes) <- storm_k562_sc_tes.common_bulk.tes$te_name_to_plot

# STORM cell number filter is needed
cell_filt <- rowSums(bulk_and_storm_sc_tes[,c(3:93)] != 0)/91
summary(cell_filt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1099  0.2527  0.3505  0.5604  1.0000 
cell_filt <- cell_filt >= 0.1
bulk_and_storm_sc_tes.cell_filt <- bulk_and_storm_sc_tes[cell_filt,]

#### KNOWN ONCOEXAPTATION WORK ####
# this requires some manual curation
# some of the genes are different as well as the initiating TE
# known oncoexap events by te-gene combos
curated_oncoexap_k562 <- read.delim("known_k562_oncoexaptation_shah_et_al_2023_natgen/curated_known_events.txt",
                                    header = TRUE)
# remove the "no" TEs
curated_oncoexap_k562.filt <- curated_oncoexap_k562[curated_oncoexap_k562$Found.in.STORM %in% "yes",]
# subset based on TE location - much easier
# will need to add in the MLT1A0-int-ENSG00000002822 - see notes
curated_oncoexap_k562.filt$Chr.TE <- gsub("chr", "", curated_oncoexap_k562.filt$Chr.TE)
curated_te_locs <- paste0(curated_oncoexap_k562.filt$Chr.TE, ":",
                          curated_oncoexap_k562.filt$Start.TE, "-",
                          curated_oncoexap_k562.filt$End.TE)
curated_te_locs <- data.frame(te_loc = curated_te_locs,
                              te_subfam = curated_oncoexap_k562.filt$Subfam)
library(GenomicRanges)
curated_te_locs.gr <- as(curated_te_locs$te_loc, "GRanges")
mcols(curated_te_locs.gr)$TE_subfam <- curated_te_locs$te_subfam

# intersect with the te-derived txps
storm_k562_sc_tes.locs <- data.frame(te_loc = paste0(storm_k562_sc_tes$TE_chr, ":",
                                                     storm_k562_sc_tes$TE_start, "-",
                                                     storm_k562_sc_tes$TE_stop),
                                     te_txp = storm_k562_sc_tes$transcript_id)
storm_k562_sc_tes.locs.gr <- as(storm_k562_sc_tes.locs$te_loc, "GRanges")
mcols(storm_k562_sc_tes.locs.gr)$TE_txp <- storm_k562_sc_tes.locs$te_txp

storm_k562_sc_tes.locs.gr.known <- subsetByOverlaps(storm_k562_sc_tes.locs.gr,
                                                    curated_te_locs.gr)
storm_k562_sc_tes.locs.gr.known_txps <- storm_k562_sc_tes.locs.gr.known$TE_txp
add_on_mlt10aint <- storm_k562_sc_tes[storm_k562_sc_tes$teprof3_gene_name %in% "MLT1A0-int-ENSG00000002822",]
storm_k562_sc_tes.locs.gr.known_txps <- c(storm_k562_sc_tes.locs.gr.known_txps,
                                          add_on_mlt10aint$transcript_id)

storm_k562_sc_tes.known_oncoexap.filt <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% storm_k562_sc_tes.locs.gr.known_txps,]

storm_k562_sc_tes.quant.known_oncoexap <- storm_k562_sc_tes.quant[storm_k562_sc_tes.quant$transcript_id %in% 
                                                                    storm_k562_sc_tes.known_oncoexap.filt$transcript_id,]

# collapse within TE
storm_k562_sc_tes.known_oncoexap.filt$TE_loc <- paste0(storm_k562_sc_tes.known_oncoexap.filt$TE_chr,
                                                       ":", storm_k562_sc_tes.known_oncoexap.filt$TE_start,
                                                       "-", storm_k562_sc_tes.known_oncoexap.filt$TE_stop)
annot_subset <- storm_k562_sc_tes.known_oncoexap.filt %>% 
  select(transcript_id, teprof3_gene_name)

storm_k562_sc_tes.quant.known_oncoexap <- storm_k562_sc_tes.quant.known_oncoexap %>%
  left_join(annot_subset, by = "transcript_id")

bulk_tes_storm_consensus_quant.known_oncoexap <- bulk_tes_storm_consensus_quant[bulk_tes_storm_consensus_quant$transcript_id %in% 
                                                                                  storm_k562_sc_tes.known_oncoexap.filt$transcript_id,]
bulk_tes_storm_consensus_quant.known_oncoexap <- bulk_tes_storm_consensus_quant.known_oncoexap %>% 
  left_join(annot_subset, by = "transcript_id")


storm_k562_sc_tes_aggregated_data.known <- storm_k562_sc_tes.quant.known_oncoexap %>%
  group_by(teprof3_gene_name, remapped) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")  # Sum over all transcript_id

bulk_tes_storm_consensus_quant.agg <- bulk_tes_storm_consensus_quant.known_oncoexap %>%
  group_by(teprof3_gene_name, sample) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")

# make the expression matrix
storm_k562_sc_tes_expression_matrix.known <- storm_k562_sc_tes_aggregated_data.known %>%
  spread(key=remapped, value=sum_tpm, fill=0)  # Convert to wide format

bulk_tes_storm_consensus_quant_expression_matrix <- bulk_tes_storm_consensus_quant.agg %>%
  spread(key=sample, value=sum_tpm, fill=0)


# need to manually sum a couple genes
# collapse to "L1PA2-ENSG00000106536" which is POU6F
pou6f.tu <- c("L1PA2-TU1391",
              "L1PA2-TU1393",
              "L1PA2-TU1397",
              "L1PA2-ENSG00000106536")
pou6f.storm <- storm_k562_sc_tes_expression_matrix.known[storm_k562_sc_tes_expression_matrix.known$teprof3_gene_name %in% 
                                                           pou6f.tu,]
pou6f.storm.sum <- c(teprof3_gene_name = "L1PA2-ENSG00000106536",
                     colSums(pou6f.storm[,c(2:97)]))
storm_k562_sc_tes_expression_matrix.known <- storm_k562_sc_tes_expression_matrix.known[!storm_k562_sc_tes_expression_matrix.known$teprof3_gene_name %in% 
                                                                                         pou6f.tu,]
storm_k562_sc_tes_expression_matrix.known <- rbind(storm_k562_sc_tes_expression_matrix.known,
                                                   pou6f.storm.sum)
pou6f.bulk <- bulk_tes_storm_consensus_quant_expression_matrix[bulk_tes_storm_consensus_quant_expression_matrix$teprof3_gene_name %in% 
                                                                 pou6f.tu,]
pou6f.bulk.sum <- c(teprof3_gene_name = "L1PA2-ENSG00000106536",
                    colSums(pou6f.bulk[,c(2:3)]))
bulk_tes_storm_consensus_quant_expression_matrix <- bulk_tes_storm_consensus_quant_expression_matrix[!bulk_tes_storm_consensus_quant_expression_matrix$teprof3_gene_name %in% 
                                                                                                       pou6f.tu,]
bulk_tes_storm_consensus_quant_expression_matrix <- rbind(bulk_tes_storm_consensus_quant_expression_matrix,
                                                          pou6f.bulk.sum)

# repeat for AluY-ENSG00000184226, AluY-TU459, AluY-TU461
# which is PCDH9
pcdh9.tu <- c("AluY-TU459",
              "AluY-TU461",
              "AluY-ENSG00000184226")
pcdh9.storm <- storm_k562_sc_tes_expression_matrix.known[storm_k562_sc_tes_expression_matrix.known$teprof3_gene_name %in% 
                                                         pcdh9.tu,]
pcdh9.storm[, 2:97] <- lapply(pcdh9.storm[, 2:97], function(x) as.numeric(as.character(x)))
pcdh9.storm.sum <- c(teprof3_gene_name = "AluY-ENSG00000184226",
                     colSums(pcdh9.storm[,c(2:97)]))
storm_k562_sc_tes_expression_matrix.known <- storm_k562_sc_tes_expression_matrix.known[!storm_k562_sc_tes_expression_matrix.known$teprof3_gene_name %in% 
                                                                                       pcdh9.tu,]
storm_k562_sc_tes_expression_matrix.known <- rbind(storm_k562_sc_tes_expression_matrix.known,
                                                   pcdh9.storm.sum)
pcdh9.bulk <- bulk_tes_storm_consensus_quant_expression_matrix[bulk_tes_storm_consensus_quant_expression_matrix$teprof3_gene_name %in% 
                                                                 pcdh9.tu,]
pcdh9.bulk[, 2:3] <- lapply(pcdh9.bulk[, 2:3], function(x) as.numeric(as.character(x)))
pcdh9.bulk.sum <- c(teprof3_gene_name = "AluY-ENSG00000184226",
                    colSums(pcdh9.bulk[,c(2:3)]))
bulk_tes_storm_consensus_quant_expression_matrix <- bulk_tes_storm_consensus_quant_expression_matrix[!bulk_tes_storm_consensus_quant_expression_matrix$teprof3_gene_name %in% 
                                                                                                       pcdh9.tu,]
bulk_tes_storm_consensus_quant_expression_matrix <- rbind(bulk_tes_storm_consensus_quant_expression_matrix,
                                                          pcdh9.bulk.sum)

storm_k562_sc_tes_expression_matrix.known[,2:97] <- lapply(storm_k562_sc_tes_expression_matrix.known[,2:97],
                                                            function(x) as.numeric(as.character(x)))
storm_k562_sc_tes_expression_matrix.known <- cbind(storm_k562_sc_tes_expression_matrix.known$teprof3_gene_name,
                                             log2(storm_k562_sc_tes_expression_matrix.known[,c(2:97)] + 1))
colnames(storm_k562_sc_tes_expression_matrix.known)[1] <- "gene_name"

bulk_tes_storm_consensus_quant_expression_matrix[,2:3] <- lapply(bulk_tes_storm_consensus_quant_expression_matrix[,2:3],
                                                           function(x) as.numeric(as.character(x)))
bulk_tes_storm_consensus_quant_expression_matrix <- cbind(bulk_tes_storm_consensus_quant_expression_matrix$teprof3_gene_name,
                                                   log2(bulk_tes_storm_consensus_quant_expression_matrix[,c(2:3)] + 1))
colnames(bulk_tes_storm_consensus_quant_expression_matrix)[1] <- "gene_name"

# add rownames
rownames(storm_k562_sc_tes_expression_matrix.known) <- storm_k562_sc_tes_expression_matrix.known$gene_name
storm_k562_sc_tes_expression_matrix.known$gene_name <- NULL

rownames(bulk_tes_storm_consensus_quant_expression_matrix) <- bulk_tes_storm_consensus_quant_expression_matrix$gene_name
bulk_tes_storm_consensus_quant_expression_matrix$gene_name <- NULL

# filter cells
storm_gene_counts <- storm_sce.filt$NumGenesExpressed
names(storm_gene_counts) <- paste0(colnames(storm_sce.filt),
                                   "_teprof3")
storm_gene_counts <- storm_gene_counts[storm_gene_counts > 9000]
storm_k562_sc_tes_expression_matrix.known <- storm_k562_sc_tes_expression_matrix.known[,colnames(storm_k562_sc_tes_expression_matrix.known) %in% c(names(storm_gene_counts),
                                                                                                                                                   "gene_name")]
storm_k562_sc_tes_expression_matrix.known$gene_name <- NULL
all_known_oncoexap <- cbind(bulk_tes_storm_consensus_quant_expression_matrix,
                            storm_k562_sc_tes_expression_matrix.known)
assay_types <- c("Bulk", "Bulk",
                 rep("STORM", 91))

library(ComplexHeatmap)
Heatmap(as.matrix(all_known_oncoexap),
        column_split = assay_types,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)                                                                                                                                 


assay_types <- c("Bulk", "Bulk",
                 rep("STORM", 91))

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

heatmap_legend_params <- list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12),
                              legend_width = unit(4, "cm"),
                              legend_height = unit(4, "cm"))

# what happens if I parse out the prop of cells expressing into quartiles?
cell_props <- rowSums(bulk_and_storm_sc_tes.cell_filt[,c(3:93)] != 0)/91

cell_prop_thresholds <- c(0.1, 0.5, 0.8, 1.0)
cell_prop_bins <- cut(cell_props,
                      breaks = cell_prop_thresholds,
                      include.lowest = TRUE,
                      right = FALSE,
                      labels = c("0.1-0.5", "0.5-0.8", "0.8-1"))

Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt),
        column_split = assay_types,
        row_split = cell_prop_bins,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

# high expression across cells
bulk_and_storm_sc_tes.cell_filt.highprop <- bulk_and_storm_sc_tes.cell_filt[cell_prop_bins %in% "0.8-1",]


Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt.highprop),
        column_split = assay_types,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

bulk_and_storm_sc_tes.cell_filt.midprop <- bulk_and_storm_sc_tes.cell_filt[cell_prop_bins %in% "0.5-0.8",]


Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt.midprop),
        column_split = assay_types,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

bulk_and_storm_sc_tes.cell_filt.lowprop <- bulk_and_storm_sc_tes.cell_filt[cell_prop_bins %in% "0.1-0.5",]


Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt.lowprop),
        column_split = assay_types,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)