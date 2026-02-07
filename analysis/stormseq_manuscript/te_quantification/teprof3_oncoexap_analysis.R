## process the TEProf3 output

## bulk total RNA

bulk_tes <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/bulk_total_sr_ont_directrna_directcdna/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")

# subset the above to the overlaps with known
# oncoexaptation events from https://pubmed.ncbi.nlm.nih.gov/36973455/
# that also overlap with a CAGE peak from Fantom5 (https://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/)
# lifted over the hg19 CAGE peaks to hg38 with UCSC liftover (v447)

bulk_tes_ovlp_known_oncoexap <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/bulk_rna_teprof3.known_oncoexaptation_k562.f5_cage_intersect.bed.gz",
                                           header = FALSE)
# NOTE: column V8 has the TE-transcript ID (e.g. TU1)
bulk_tes.filt <- bulk_tes[bulk_tes$transcript_id %in% bulk_tes_ovlp_known_oncoexap$V8,]

# one strategy described in the TEProf3 docs is to grab things
# that have non-zero perfect_SJ_uniqlymapped_read_downstream in the
# teprof3_output_quantification.TE.tsv.gz 
# To determine whether a TE-derived transcript is present or not in one sample, 
# we recommned to use the following three criteria: 
# (1) > 1TPM 
# (2) >= 1 perfect_SJ_uniqlymapped_read_downstream 
# (3) SJ_uniqlymapped_read_upstream <= 0.5 * perfect_SJ_uniqlymapped_read_downstream
# NOTE: Josh mentioned upping the TPM threshold significantly to find the most highly expressed events
bulk_tes_quant <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/bulk_total_sr_ont_directrna_directcdna/teprof3_output_quantification.TE.tsv.gz")

# filter based on recommended TEProf3 docs
bulk_tes_quant.filt <- bulk_tes_quant[bulk_tes_quant$sum_tpm >= 1,]
bulk_tes_quant.filt <- bulk_tes_quant.filt[bulk_tes_quant.filt$perfect_SJ_uniqlymapped_read_downstream > 0,]
bulk_tes_quant.filt <- bulk_tes_quant.filt[bulk_tes_quant.filt$SJ_uniqlymapped_read_upstream <= (0.5 * bulk_tes_quant.filt$perfect_SJ_uniqlymapped_read_downstream),]
table(length(unique(bulk_tes_quant.filt$transcript_id)))
# 2197

# pull in the bulk total RNA short-read data that has been run in guided mode
# using the STORM K-562 pseudobulk consensus TE GTF
bulk_tes_storm_consensus <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/bulk_with_storm_consensus_gtf/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")

bulk_tes_storm_consensus_quant <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/bulk_with_storm_consensus_gtf/teprof3_output_quantification.TE.tsv.gz")

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

plot(hist(bulk_tes_storm_consensus_quant.filt$sum_tpm))

# tpm_filt <- bulk_tes_storm_consensus_quant.filt$sum_tpm > 50
# bulk_tes_storm_consensus_quant.filt.tpm_filt <- bulk_tes_storm_consensus_quant.filt[tpm_filt,]
# table(length(unique(bulk_tes_storm_consensus_quant.filt.tpm_filt$transcript_id)))
# 100

## STORM
storm_k562_sc_tes <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/storm_deconv/teprof3_output_filter_transcript_TE_transcript_consensus.tsv")
storm_k562_sc_tes.quant <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/teprof3/storm_deconv/teprof3_output_quantification.TE.tsv.gz")

# pull in the expression SCE from the benchmarking to filter to QC'd cells
storm_kb_res <- readRDS("~/Documents/manuscripts/storm_seq/subsample_analysis/complexity_curves/storm/storm_raw_complexity_curve_gene_txp_counts_ens101_with_erccs.rds")
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
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

# stick to 1M subsample
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

# need to remap the storm te column names....
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
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
  filter(nonzero_cells / total_cells >= 0.1) %>%  # Keep if at least 20% have non-zero sum_tpm
  # filter(perfect_sj_nonzero_cells / total_cells >= 0.1) %>% # Think of this like a second layer of evidence if TPM is low
  ungroup()
length(unique(storm_k562_sc_tes.quant.filt$transcript_id))
# 663

# relaxed
# this is for looking at associations with fusions
storm_k562_sc_tes.quant.filt.relaxed <- storm_k562_sc_tes.quant %>%
  group_by(transcript_id) %>%
  mutate(total_cells = n(),                      # Count total occurrences of each transcript_id
         nonzero_cells = sum(stringtie_tpm > 0),   # Count non-zero stringtie_tpm values
         median_tpm = median(stringtie_tpm)) %>%     # Calculate the median stringtie_tpm per transcript
  mutate(perfect_sj_nonzero_cells = sum(perfect_SJ_uniqlymapped_read_downstream > 0)) %>%
  filter(nonzero_cells / total_cells >= 0.2 &
           nonzero_cells / total_cells <= 0.9) %>%   # Keep transcripts with 20%-80% non-zero expression
  # filter(median_tpm >= 0.1) %>%                         # Filter out transcripts with median TPM <= 1 (adjust as needed)
  ungroup()

# how many TE txps
# length(unique(storm_k562_sc_tes.quant.filt.relaxed$transcript_id))
# 327 TE-derived transcripts

# how many unique genes
# length(unique(storm_k562_sc_tes.quant.filt.relaxed$gene_name))
# 303 unique genes

# are these found in the bulk?
# storm_single_cell_oncoexaptations <- unique(storm_k562_sc_tes.quant.filt$transcript_id)
# table(storm_single_cell_oncoexaptations %in% bulk_tes_storm_consensus_quant.filt$transcript_id)
# FALSE  TRUE 
# 25   437

# let's just subset to bulk TE-txps since that's our "ground truth"
storm_single_cell_oncoexaptations.bulk <- storm_k562_sc_tes.quant[storm_k562_sc_tes.quant$transcript_id %in% bulk_tes_storm_consensus_quant.filt$transcript_id,]

# fusion associations
# storm_single_cell_oncoexaptations.bulk <- storm_k562_sc_tes.quant.filt.relaxed[storm_k562_sc_tes.quant.filt.relaxed$transcript_id %in% bulk_tes_storm_consensus_quant.filt$transcript_id,]
# storm_k562_sc_tes.quant.filt.relaxed <- storm_single_cell_oncoexaptations.bulk

# how many genes?
# length(unique(storm_k562_sc_tes.quant.filt$gene_name))
# 229

# what do the non-bulk ones look like?
# non_bulk <- storm_single_cell_oncoexaptations[!storm_single_cell_oncoexaptations %in% bulk_tes_storm_consensus_quant.filt$transcript_id]
# storm_k562_sc_tes.quant.filt.non_bulk <- storm_k562_sc_tes.quant.filt[storm_k562_sc_tes.quant.filt$transcript_id %in% non_bulk,]

# are there any that are expressed in _all_ of the cells?
# storm_k562_sc_tes.quant.filt_all <- storm_k562_sc_tes.quant %>%
#   group_by(transcript_id) %>%
#   mutate(total_cells = n(),  # Count total occurrences of each transcript_id
#          nonzero_cells = sum(sum_tpm > 0)) %>%  # Count non-zero `sum_tpm` values
#   mutate(perfect_sj_nonzero_cells = sum(perfect_SJ_uniqlymapped_read_downstream > 0)) %>%
#   filter(nonzero_cells / total_cells >= 0.9) %>%  # Keep if at least 90% have non-zero sum_tpm
#   filter(perfect_sj_nonzero_cells / total_cells >= 0.1) %>%
#   ungroup()
# length(unique(storm_k562_sc_tes.quant.filt_all$transcript_id))
# 140!!

# are these found in the bulk?
# storm_all_cell_oncoexaptations <- unique(storm_k562_sc_tes.quant.filt_all$transcript_id)
# table(storm_all_cell_oncoexaptations %in% bulk_tes_storm_consensus_quant.filt$transcript_id)
# TRUE 
# 140 

# what are the genes? 34 of them
# unique(storm_k562_sc_tes.quant.filt_all$gene_name)
# [1] "ENSG00000287756" "ENSG00000162594" "TU457"           "ENSG00000069702" "ENSG00000223756"
# [6] "ENSG00000167355" "ENSG00000251381" "ENSG00000174718" "ENSG00000111300" "ENSG00000257636"
# [11] "ENSG00000196547" "ENSG00000275016" "ENSG00000260289" "ENSG00000262879" "ENSG00000230258"
# [16] "ENSG00000132204" "ENSG00000263711" "ENSG00000265843" "ENSG00000079999" "ENSG00000266976"
# [21] "ENSG00000267243" "ENSG00000286769" "ENSG00000235335" "ENSG00000233922" "TU12272"        
# [26] "ENSG00000227706" "ENSG00000237643" "ENSG00000106397" "ENSG00000230257" "ENSG00000161040"
# [31] "ENSG00000228742" "ENSG00000064419" "ENSG00000105939" "ENSG00000234722" "ENSG00000253981"
# [36] "ENSG00000035681" "ENSG00000253394" "ENSG00000229140" "ENSG00000095261" "ENSG00000237311"
# [41] "ENSG00000181433"

# filter the annotations to learn a bit more about these
# storm_k562_sc_tes.filt_all <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% storm_all_cell_oncoexaptations,]
# table(storm_k562_sc_tes.filt_all$TE_family)
# Alu        ERV1        ERVL   ERVL-MaLR hAT-Charlie          L1          L2         MIR 
# 26          57           8          35           1           6           4           3 

# alright... let's plot them
library(ComplexHeatmap)
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

# fusion associated hunting
table(mono_exons$transcript_id %in% storm_k562_sc_tes.quant.filt.relaxed$transcript_id)
# FALSE 
# 12930 

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

# add a TPM filter for the bulk
# bulk_k562_sc_tes_all_expression_matrix.log.tpm_filt <- bulk_k562_sc_tes_all_expression_matrix.log[rowMeans(bulk_k562_sc_tes_all_expression_matrix.log[,c(2:3)]) >= 2,]

#bulk_and_storm_sc_tes <- inner_join(bulk_k562_tes_expression_matrix, storm_k562_sc_tes_expression_matrix, by = "gene_name")
bulk_and_storm_sc_tes <- inner_join(bulk_k562_tes_expression_matrix,
                                    storm_k562_sc_tes_expression_matrix, by = "gene_name")
rownames(bulk_and_storm_sc_tes) <- bulk_and_storm_sc_tes$gene_name
bulk_and_storm_sc_tes$gene_name <- NULL
# bulk_and_storm_sc_tes$storm_pseudobulk <- rowMeans(bulk_and_storm_sc_tes[,c(4:99)])

# go through and find the TE-Gene combo for annotations
# e.g. L2c-ENSG00000229140
# pull in the GTF and find out what genes we have going on
gtf <- rtracklayer::import("~/Documents/manuscripts/storm_seq/rrna/Homo_sapiens.GRCh38.101.ercc92patched.gtf.gz")
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

# ENSG00000229140
# ens_4 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000229140",]
# ens_4 <- ens_4[ens_4$percentage_of_expression > 0,]
# table(ens_4$transcript_id)
# # TU15003 TU15035 
# # 66      27 
# tu_check <- unique(ens_4$transcript_id)
# ens_4 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$transcript_id %in% tu_check,]
# tu_check.annot <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% tu_check,]
# ens_4 <- ens_4[ens_4$percentage_of_expression > 50,]
# table(ens_4$transcript_id)
# # TU15035 is unique!
# tu_keep <- c(tu_keep, "TU15035")
# # Upon further inspection, TU15003 is it's own thing!
# tu_keep <- c(tu_keep, "TU15003")

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

# ENSG00000234948
# ens_6 <- storm_k562_sc_tes.common_bulk_dedupe.fine[storm_k562_sc_tes.common_bulk_dedupe.fine$gene_name %in% "ENSG00000234948",]
# ens_6 <- ens_6[ens_6$percentage_of_expression > 0,]
# table(ens_6$transcript_id)
# # TU8327 TU8328
# # 23      6
# 
# # TU8327 is dominant
# tu_keep <- c(tu_keep, "TU8327")

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
curated_oncoexap_k562 <- read.delim("~/Documents/manuscripts/storm_seq/te_quant/known_k562_oncoexaptation_shah_et_al_2023_natgen/curated_known_events.txt",
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


# convert the gene names# convert the gene names
# known_gene_names <- storm_k562_sc_tes_expression_matrix.known$gene_name
# gtf.genes.known_oncoexap.m <- match(known_gene_names,
#                                     gtf.genes.known_oncoexap$gene_id)
# gtf.genes.known_oncoexap <- gtf.genes.known_oncoexap[gtf.genes.known_oncoexap.m,]
# known_oncoexap_k562.m <- match(gtf.genes.known_oncoexap$gene_name,
#                                known_oncoexap_k562$V2)
# known_oncoexap_k562 <- known_oncoexap_k562[known_oncoexap_k562.m,]
# rownames(storm_k562_sc_tes_expression_matrix.known) <- paste0(known_oncoexap_k562$V1,
#                                                               "-",
#                                                               known_oncoexap_k562$V2)
# storm_k562_sc_tes_expression_matrix.known$gene_name <- NULL
# 
# # add bulk data now for comparison
# bulk_tes_storm_consensus_quant_expression_matrix <- bulk_tes_storm_consensus_quant_expression_matrix[bulk_tes_storm_consensus_quant_expression_matrix$gene_name %in% gtf.genes.known_oncoexap$gene_id,]
# rownames(bulk_tes_storm_consensus_quant_expression_matrix) <- paste0(known_oncoexap_k562$V1,
#                                                                      "-",
#                                                                      known_oncoexap_k562$V2)
# bulk_tes_storm_consensus_quant_expression_matrix$gene_name <- NULL
# plot(rowMeans(bulk_tes_storm_consensus_quant_expression_matrix),
#      rowMeans(storm_k562_sc_tes_expression_matrix.known))

# Exclude any cell with less than 9k genes
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

Heatmap(as.matrix(all_known_oncoexap),
        column_split = assay_types,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)                                                                                                                                 


# table(known_oncoexap_k562.tegene %in% rownames(bulk_and_storm_sc_tes)) 

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
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1099  0.2170  0.3407  0.4385  0.6484  1.0000 
# cell_props.quartiles <- cut(cell_props, 
#                             breaks = quantile(cell_props,
#                                               probs = seq(0, 1, 0.25),
#                                               na.rm = TRUE), 
#                             include.lowest = TRUE, 
#                             labels = c("Q1", "Q2", "Q3", "Q4"))

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

# can I overlay bcr-abl and nup214-xkr3 here?
# go down and use the binary fusion matrix for labeling
bcr_abl_fusions <- binary_mat_fusions$`BCR--ABL1`
names(bcr_abl_fusions) <- binary_mat_fusions$remapped
bcr_abl_fusions <- bcr_abl_fusions[names(bcr_abl_fusions) %in% colnames(bulk_and_storm_sc_tes.cell_filt)]

nup_xkr_fusions <- binary_mat_fusions$`NUP214--XKR3`
names(nup_xkr_fusions) <- binary_mat_fusions$remapped
nup_xkr_fusions <- nup_xkr_fusions[names(nup_xkr_fusions) %in% colnames(bulk_and_storm_sc_tes.cell_filt)]

bulk_and_storm_sc_tes.cell_filt.fusion_filt <- bulk_and_storm_sc_tes.cell_filt[,c(TRUE, TRUE, colnames(bulk_and_storm_sc_tes.cell_filt)[3:93]
                                                                               %in% names(bcr_abl_fusions))]
# add in the fusions
bcr_abl_fusions.bulk <- c("bulk_total_k562_rep1" = 1,
                          "bulk_total_k562_rep2" = 1)
bcr_abl_fusion.all <- c(bcr_abl_fusions.bulk,
                        bcr_abl_fusions)
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

# pull out LTR1A2-PURPL
purpl <- bulk_and_storm_sc_tes.cell_filt.midprop[grep("PURPL", rownames(bulk_and_storm_sc_tes.cell_filt.midprop)),]
purpl_cells <- purpl[,c(3:93)]
purpl_cells <- purpl_cells[,purpl_cells > 0]

# how many cells have bcr-abl?
bcrabl <- binary_mat_fusions[,c(1:2)]
bcrabl.cells <- bcrabl$remapped[bcrabl$`BCR--ABL1` > 0]
table(colnames(purpl_cells) %in% bcrabl.cells)
# FALSE  TRUE 
# 20    26 
# how many cells have nup214-xkr3?
nup <- binary_mat_fusions[,c(1,3)]
nup.cells <- nup$remapped[nup$`NUP214--XKR3` > 0]
table(colnames(purpl_cells) %in% nup.cells)
# FALSE  TRUE 
# 24    22
# split!
# meh - move on
# grab a couple purpl expressing cells to plot
# K15 and A17
# original K3 and K9 look good

bulk_and_storm_sc_tes.cell_filt.lowprop <- bulk_and_storm_sc_tes.cell_filt[cell_prop_bins %in% "0.1-0.5",]


Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt.lowprop),
        column_split = assay_types,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

# Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt),
#         column_split = assay_types,
#         show_row_names = FALSE,
#         show_column_names = FALSE,
#         col = getJetColors(circlize.cols = TRUE),
#         name = "log2(TPM+1)",
#         heatmap_legend_param = heatmap_legend_params)

plot(hist(rowMeans(bulk_and_storm_sc_tes.cell_filt[,c(1:2)])))

discovery_te_oncoexap <- data.frame(te_oncoexaptation = rownames(bulk_and_storm_sc_tes.cell_filt),
                                    bulk_mean_exp = rowMeans(bulk_and_storm_sc_tes.cell_filt[,c(1:2)]),
                                    storm_mean_exp = rowMeans(bulk_and_storm_sc_tes.cell_filt[,c(3:93)]),
                                    proportion_storm_cells_expressed = rowSums(bulk_and_storm_sc_tes.cell_filt[,c(3:93)] != 0)/91)

# add in a TPM filter for high expression bulk TE-oncoexaptation events
bulk_and_storm_sc_tes.cell_filt.tpm_filt <- bulk_and_storm_sc_tes.cell_filt[rowMeans(bulk_and_storm_sc_tes.cell_filt[,c(1:2)]) >= 3,]

# add cell filtration filter again
cell_filt <- rowSums(bulk_and_storm_sc_tes.cell_filt.tpm_filt[,c(3:93)] != 0)/91
summary(cell_filt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2198  0.3736  0.7253  0.6504  0.9011  1.0000 
cell_filt <- cell_filt >= 0.5
bulk_and_storm_sc_tes.cell_filt.tpm_filt <- bulk_and_storm_sc_tes.cell_filt.tpm_filt[cell_filt,]

# filter out TU and antisense txps (e.g., AS1)
# this was unstranded so hard to know which strand it came from
# as the anti-sense txps often overlap the sense version
bulk_and_storm_sc_tes.cell_filt.tpm_filt <- bulk_and_storm_sc_tes.cell_filt.tpm_filt[grep("-TU|-AS1", rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt), invert = TRUE),]

# see if we can remap some of the AC and AL genes...
txp_ids <- rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)
names(txp_ids) <- gsub(".*-", "", txp_ids)

gene_ids <- gtf.genes[gtf.genes$gene_name %in% names(txp_ids),]$gene_id
names(gene_ids) <- gtf.genes[gtf.genes$gene_name %in% names(txp_ids),]$gene_name

gene_ids.m <- match(names(txp_ids),
                    names(gene_ids))
gene_ids <- gene_ids[gene_ids.m]

# pull out the AC, AF, AP, AL genes
gene_ids.sub <- gene_ids[grep("^AC|^AF|^AP|^AL", names(gene_ids))]
gene_ids.sub <- gene_ids.sub[grep("\\.1$|\\.2$|\\.3$|\\.5$", names(gene_ids.sub))]
gene_ids.sub
# AC097528.1        AL713998.1        AL365226.1        AL138916.1 
# "Lnc-MFSD8-6"    "Lnc-LMBRD1-10"    "EVADR"           "Lnc-ADGB-1" 
# AL034397.2        AC079949.5        AC015574.1        AC109446.3 
# "Lnc-HEPH-2"      "Lnc-TMEM132C-6" "Lnc-NR2F2-14/16" "XYLT1-AS" 
# AC005208.1        AP001531.1        AC079466.1        AC005381.1 
# "Lnc-ABCA5-2"     "HSALNG0145415"   "ERVE-5"         "Lnc-UQCRFS1-3" 
# AC110614.1        AF127936.1          AP000561.1 
# "SUCLA2P2"       "Lnc-SAMSN1-1/2/3"   "Lnc-NCAM2-4/13" 


# not too many to fix... just do it by hand
# AC097528.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC097528.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MER52A-Lnc-MFSD8-6"
# AL713998.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AL713998.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "LTR7C-Lnc-LMBRD1-10"
# AL365226.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AL365226.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MER48-EVADR"
# AL138916.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AL138916.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MLT2A2-Lnc-ADGB-1"
# AL034397.2
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AL034397.2",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MIR3-Lnc-HEPH-2"
# AC079949.5
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC079949.5",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "LTR1A2-Lnc-TMEM132C-6"
# AC015574.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC015574.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "AluYm1-Lnc-NR2F2-14/16"
# AC109446.3
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC109446.3",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "AluY-XYLT1-AS"
# AC005208.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC005208.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MER48-Lnc-ABCA5-2"
# AP001531.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AP001531.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "LTR2-HSALNG0145415"
# AC079466.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC079466.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "LTR2C-ERVE-5"
# AC005381.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC005381.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "MLT1K-Lnc-UQCRFS1-3"
# AC110614.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AC110614.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "LTR17-SUCLA2P2"
# AF127936.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AF127936.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "L1MB1-Lnc-SAMSN1-1/2/3"
# AP000561.1
rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt)[grep("AP000561.1",
                                                        rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt))] <- "AluY-Lnc-NCAM2-4/13"

# now get rid of the rest of the AC/AL/AP/AF txps as we just don't know
# what they are and this will not add to the final figure


Heatmap(as.matrix(bulk_and_storm_sc_tes.cell_filt.tpm_filt),
        column_split = assay_types,
        show_row_names = TRUE,
        show_column_names = FALSE,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)


curated_te_oncoexap <- data.frame(te_oncoexaptation = rownames(bulk_and_storm_sc_tes.cell_filt.tpm_filt),
                                  bulk_mean_exp = rowMeans(bulk_and_storm_sc_tes.cell_filt.tpm_filt[,c(1:2)]),
                                  storm_mean_exp = rowMeans(bulk_and_storm_sc_tes.cell_filt.tpm_filt[,c(3:93)]),
                                  proportion_storm_cells_expressed = rowSums(bulk_and_storm_sc_tes.cell_filt.tpm_filt[,c(3:93)] != 0)/91)




# NOTES:
# Interesting TE-derived transcripts
# lncRNA CCDC26 has 4 separate TE-derived transcripts that are all the same...
# likely due to the assembly/stitching process. 
# TU15030 TU15032 TU15034 TU15035 
# 37      52      49      27 
# all have the same start with minor variations in ends
# the TE-txp is also the same: L2c-ENSG00000229140

# what is the pseudobulk expression correlation with mean of bulk
bulk_exp_mean <- rowMeans(bulk_and_storm_sc_tes[,c(1:2)])
storm_exp_mean <- rowMeans(bulk_and_storm_sc_tes[,c(3:93)])

aves <- data.frame(bulk = bulk_exp_mean,
                   storm = storm_exp_mean,
                   genes = rownames(bulk_and_storm_sc_tes))
library(ggplot2)
ggplot(aves, aes(x = bulk,
                 y = storm)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlim(c(0, 10)) +
  ylim(c(0, 10))

cor(bulk_exp_mean,
    storm_exp_mean)
# 0.6720266

# pull in the GTF and find out what genes we have going on
gtf <- rtracklayer::import("~/Documents/manuscripts/storm_seq/rrna/Homo_sapiens.GRCh38.101.ercc92patched.gtf.gz")
gtf.genes <- gtf[mcols(gtf)$type == "gene",]
gtf.genes.tes <- gtf.genes[gtf.genes$gene_id %in% unique(storm_k562_sc_tes.quant.filt_all$gene_name),]
gtf.genes.tes$gene_name
# [1] "AL031293.1"  "IL23R"       "TGFBR3"      "AC064869.1"  "B3GALT1-AS1" "AL713998.1"  "AL365226.1" 
# [8] "PLOD3"       "FBXL13"      "NFE4"        "LINC02577"   "TNPO3"       "ZC3HAV1"     "LINC01287"  
# [15] "AL034397.2"  "SAGE1"       "ALG1L13P"    "NSMAF"       "LINC00534"   "CCDC26"      "PSMD5"      
# [22] "TSSC2"       "OR51B5"      "LINC00958"   "RESF1"       "NAA25"       "G2E3-AS1"    "MAN2A2"     
# [29] "AC015574.1"  "AC093515.1"  "AC005670.3"  "AC005208.1"  "LINC00470"   "LINC02864"   "LINC01029"  
# [36] "KEAP1"       "AC005381.1"  "AC079466.1"  "LINC01694"   

# look in the STORM deconv now after filtering
gtf.genes.storm_deconv_tes <- gtf.genes[gtf.genes$gene_id %in% unique(storm_k562_sc_tes.quant.filt$gene_name),]
gtf.genes.storm_deconv_tes$gene_name
# [1] "AL732372.2"   "MTHFR"        "KAZN-AS1"     "AL031293.1"   "IL23R"        "AL359894.1"  
# [7] "TGFBR3"       "CCDC18"       "CDC14A"       "SLC16A1-AS1"  "LINC02805"    "RRNAD1"      
# [13] "DCAF6"        "BPNT1"        "ZNF670"       "ZNF124"       "AC108488.1"   "RAB10"       
# [19] "AC006369.1"   "FAM178B"      "TGFBRAP1"     "LINC01885"    "AC064869.1"   "B3GALT1-AS1" 
# [25] "CCDC173"      "AC010894.2"   "AC096555.1"   "LINC01090"    "PPIL3"        "PTH2R"       
# [31] "AC090004.2"   "ZNF385D-AS1"  "C3orf86"      "MRPS22"       "SEC62"        "AC073365.1"  
# [37] "AFAP1-AS1"    "AC097512.1"   "LINC02506"    "LINC02616"    "CXCL11"       "ARHGAP24"    
# [43] "DAPP1"        "LINC02432"    "MIR3945HG"    "AC106771.1"   "LINC01194"    "GTF2H2B"     
# [49] "LIX1-AS1"     "LINC02163"    "LINC01184"    "AC112178.1"   "AC010378.1"   "SPDL1"       
# [55] "NQO2"         "AL021978.1"   "CMTR1"        "AL713998.1"   "AL365226.1"   "CASP8AP2"    
# [61] "PKIB"         "RNF217"       "LINC02528"    "AL138916.1"   "AL033504.1"   "LINC02840"   
# [67] "NUDT1"        "AC007128.2"   "SCIN"         "PLOD3"        "FBXL13"       "NFE4"        
# [73] "LINC02577"    "IMMP2L"       "LINC02476"    "TNPO3"        "ZC3HAV1"      "LINC01287"   
# [79] "AL034397.2"   "KCTD9P2"      "TENM1"        "FIRRE"        "CT45A10"      "SAGE1"       
# [85] "FGF13"        "GABRE"        "AC246817.2"   "AC246817.1"   "ALG1L13P"     "AC100861.1"  
# [91] "NSMAF"        "LINC00534"    "ANGPT1"       "LINC01608"    "LINC00861"    "CCDC26"      
# [97] "MINCR"        "DOCK8"        "GALNT12"      "PSMD5"        "ENTR1"        "PMPCA"       
# [103] "TSSC2"        "HBE1"         "OR51B5"       "LINC00958"    "LTO1"         "MRPL48"      
# [109] "TAF1D"        "DCLRE1C"      "VIM-AS1"      "ST8SIA6-AS1"  "SELENOOLP"    "ZNF239"      
# [115] "TMEM273"      "ADAMTS14"     "EXOC6"        "SOX5"         "RESF1"        "AC025575.2"  
# [121] "CRADD"        "NT5DC3"       "NAA25"        "LINC01234"    "CFAP251"      "AC156455.1"  
# [127] "AC079949.5"   "TMEM132D-AS1" "EXOSC8"       "UTP14C"       "PCDH9"        "LINC00383"   
# [133] "AL445255.1"   "LINC02327"    "G2E3-AS1"     "LINC00520"    "AC016526.2"   "ESRRB"       
# [139] "LINC00239"    "AC100834.2"   "TMEM62"       "DMXL2"        "PIGB"         "LINC01169"   
# [145] "MAN2C1"       "MAN2A2"       "AC015574.1"   "WASH3P"       "AC093515.1"   "HNRNPCP4"    
# [151] "AC109446.3"   "ACSM3"        "PSME3IP1"     "LINC02193"    "NLRP1"        "UBBP4"       
# [157] "MYO19"        "AC005670.3"   "TRIM25"       "POLG2"        "AC005208.1"   "CCDC57"      
# [163] "LINC00470"    "AP001531.1"   "LINC02864"    "LINC01029"    "FASTKD5"      "LINC01524"   
# [169] "AC009005.1"   "ZNF317"       "KEAP1"        "ZNF44"        "UCA1"         "AC005381.1"  
# [175] "AC079466.1"   "ERCC1"        "CA11"         "FUT1"         "ZNF83"        "ZNF550"      
# [181] "AF127936.1"   "MIR548XHG"    "LINC01687"    "AP000561.1"   "DSCR4"        "DSCR8"       
# [187] "LINC01694" 

gtf.genes.storm_deconv_tes[seqnames(gtf.genes.storm_deconv_tes) %in% "22",]

## see if the TU15566 (found in ABL1) are found in translocated cells or not
storm_gene_fusions <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/storm/storm_star_fusion_full_depth_k562_agg_deconv.fusions.abridged.tsv")
storm_gene_fusions.bcrabl <- storm_gene_fusions[storm_gene_fusions$X.FusionName %in% "BCR--ABL1",]
storm_gene_fusions.nup214xkr3 <- storm_gene_fusions[storm_gene_fusions$X.FusionName %in% "NUP214--XKR3",]

# going to need to remap these....
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
                          header = FALSE)

# remap in a vectorized fashion
storm_remap <- setNames(storm_remap$V2, storm_remap$V1)

# BCR--ABL1
storm_gene_fusions.bcrabl$remapped <- storm_remap[storm_gene_fusions.bcrabl$Cell]
storm_gene_fusions.bcrabl.cells <- unique(storm_gene_fusions.bcrabl$remapped)

storm_TU15566 <- storm_k562_sc_tes.quant.filt.relaxed[storm_k562_sc_tes.quant.filt.relaxed$transcript_id %in% "TU15566",]
storm_TU15566 <- storm_TU15566[storm_TU15566$stringtie_tpm > 0,]
table(gsub("_teprof3", "", storm_TU15566$remapped) %in% storm_gene_fusions.bcrabl.cells)
# FALSE  TRUE 
# 3     7 


# NUP214--XKR3
storm_gene_fusions.nup214xkr3$remapped <- storm_remap[storm_gene_fusions.nup214xkr3$Cell]
storm_gene_fusions.nup214xkr3.cells <- unique(storm_gene_fusions.nup214xkr3$remapped)

storm_TU8707 <- storm_k562_sc_tes.quant.filt.relaxed[storm_k562_sc_tes.quant.filt.relaxed$transcript_id %in% "TU8707",]
storm_TU8707 <- storm_TU8707[storm_TU8707$stringtie_tpm > 0,]
table(gsub("_teprof3", "", storm_TU8707$remapped) %in% storm_gene_fusions.nup214xkr3.cells)

# look in all
storm_gene_fusions$remapped <- storm_remap[storm_gene_fusions$Cell]
storm_gene_fuxions.nonbcrabl <- c("F3", "D8", "A14")
storm_gene_fusions.sub <- storm_gene_fusions[storm_gene_fusions$remapped %in% storm_gene_fuxions.nonbcrabl,]

## run chi-square for associations with fusions
fusions_of_interest <- c("BCR--ABL1", "NUP214--XKR3",
                         "BAG6--SLC44A4", "XACT--LRCH2",
                         "C16orf87--ORC6", "UPF3A--CDC16",
                         "IMMP2L--DOCK4")
fusion_partners <- c("BCR", "ABL1",
                     "NUP214", "XKR3",
                     "BAG6", "SLC44A4",
                     "XACT", "LRCH2",
                     "C16orf87", "ORC6",
                     "UPF3A", "CDC16",
                     "IMMP2L", "DOCK4")
# gather the Ensembl Gene IDs for all of these
fusion_partners.gtf <- gtf.genes[gtf.genes$gene_name %in% fusion_partners,]
fusion_partners_ens_ids <- fusion_partners.gtf$gene_id
names(fusion_partners_ens_ids) <- fusion_partners.gtf$gene_name

# gather all TE-derived transcripts with these genes in them
# NOTE: the 'gene_name' column is incomplete and need to use 'tss_gene_name'
# where the implication is that the TSS of the TE-derived transcript
# is found within that gene...
storm_te_txps_by_fusion_partner <- storm_k562_sc_tes[storm_k562_sc_tes$tss_gene_name %in% fusion_partners_ens_ids,]
storm_te_txps_by_fusion_partner.txp_ids <- storm_te_txps_by_fusion_partner$transcript_id
storm_te_txps_by_fusion_partner.txp_ids <- setNames(storm_te_txps_by_fusion_partner$tss_gene_name,
                                                    storm_te_txps_by_fusion_partner$transcript_id)

storm_te_txps_by_fusion_partner.quant <- storm_k562_sc_tes.quant[storm_k562_sc_tes.quant$transcript_id %in% storm_te_txps_by_fusion_partner$transcript_id,]
storm_te_txps_by_fusion_partner.quant$gene_name <- storm_te_txps_by_fusion_partner.txp_ids[storm_te_txps_by_fusion_partner.quant$transcript_id]

# aggregate stringtie_tpm per gene_name and sample
storm_te_txps_by_fusion_partner_agg <- storm_te_txps_by_fusion_partner.quant %>%
  group_by(gene_name, remapped) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")  # Sum over all transcript_id

# make the expression matrix
storm_te_txps_by_fusion_partner_exp_matrix <- storm_te_txps_by_fusion_partner_agg %>%
  spread(key=remapped, value=sum_tpm, fill=0)
storm_te_txps_by_fusion_partner_exp_matrix <- as.data.frame(storm_te_txps_by_fusion_partner_exp_matrix)

# make it binary for chi-sq test# make it binary for chi-sq test
rownames(storm_te_txps_by_fusion_partner_exp_matrix) <- storm_te_txps_by_fusion_partner_exp_matrix$gene_name
storm_te_txps_by_fusion_partner_exp_matrix$gene_name <- NULL
storm_te_txps_by_fusion_partner_exp_matrix.binary <- storm_te_txps_by_fusion_partner_exp_matrix
storm_te_txps_by_fusion_partner_exp_matrix.binary[storm_te_txps_by_fusion_partner_exp_matrix.binary > 0] <- 1

# make the gene names easier to query
rownames(storm_te_txps_by_fusion_partner_exp_matrix.binary)
# [1] "ENSG00000097007" "ENSG00000128512" "ENSG00000172967"
# [4] "ENSG00000184903" "ENSG00000241743"

rownames(storm_te_txps_by_fusion_partner_exp_matrix.binary) <- c("ABL1",
                                                                 "DOCK4",
                                                                 "XKR3",
                                                                 "IMMP2L",
                                                                 "XACT")


# aggregate fusions per cell
# Filter only cells that contain at least one fusion of interest
filtered_fusions <- storm_gene_fusions %>%
  filter(X.FusionName %in% fusions_of_interest)

# Create binary matrix
binary_mat_fusions <- filtered_fusions %>%
  mutate(Present = 1) %>%
  group_by(remapped, X.FusionName) %>%  # Ensure one row per `remapped` + `X.FusionName`
  summarise(Present = max(Present), .groups = "drop") %>%  # Collapse duplicates
  pivot_wider(names_from = X.FusionName, values_from = Present, values_fill = 0)

# Ensure all fusions_of_interest are columns (if some are missing)
binary_mat_fusions <- binary_mat_fusions %>%
  select(remapped, all_of(fusions_of_interest)) %>%
  replace(is.na(.), 0)
binary_mat_fusions$remapped <- paste0(binary_mat_fusions$remapped,
                                      "_teprof3")
table(colnames(storm_te_txps_by_fusion_partner_exp_matrix.binary) %in%
        binary_mat_fusions$remapped)
# FALSE  TRUE 
# 9    87

# subset to avoid missing data
storm_te_txps_by_fusion_partner_exp_matrix.binary <- storm_te_txps_by_fusion_partner_exp_matrix.binary[,colnames(storm_te_txps_by_fusion_partner_exp_matrix.binary) %in% binary_mat_fusions$remapped]
all(colnames(storm_te_txps_by_fusion_partner_exp_matrix.binary) %in% binary_mat_fusions$remapped)
# TRUE

# chi-square tests
# Get the list of gene names (from gene expression matrix)
gene_list <- rownames(storm_te_txps_by_fusion_partner_exp_matrix.binary)  # e.g., "ABL1", "DOCK4", "XKR3", "IMMP2L", "XACT"

# Extract both fusion partners from fusion column names
fusion_partners <- colnames(binary_mat_fusions)[-1] %>%
  strsplit("--") %>%
  lapply(function(x) c(x[1], x[2])) %>% 
  setNames(colnames(binary_mat_fusions)[-1])

# Store results
fisher_results <- list()
p_values <- c()
odds_ratios <- c()

# Perform chi-square test for each gene in gene_list
for (gene in gene_list) {
  
  # Identify fusions where the gene appears in either position
  related_fusions <- names(fusion_partners)[sapply(fusion_partners, function(x) gene %in% x)]
  
  if (length(related_fusions) == 0) {
    next  # Skip if no fusions involve this gene
  }
  
  # Extract presence/absence of fusions for this gene
  fusion_status <- binary_mat_fusions %>%
    select(remapped, all_of(related_fusions)) %>%
    mutate(Fusion_Present = rowSums(across(-remapped)) > 0) %>%  # Any fusion involving this gene
    select(remapped, Fusion_Present)
  
  # Directly extract pre-binarized gene expression values
  gene_status <- data.frame(
    remapped = colnames(storm_te_txps_by_fusion_partner_exp_matrix.binary),
    Gene_Expressed = as.numeric(storm_te_txps_by_fusion_partner_exp_matrix.binary[gene, ])  # Already binary (0/1)
  )
  
  # Merge data
  merged_data <- inner_join(fusion_status, gene_status, by = "remapped")
  
  # Create contingency table
  contingency_table <- table(merged_data$Fusion_Present, merged_data$Gene_Expressed)
  print(contingency_table)
  
  # Perform Fisher’s Exact Test
  fisher_test <- fisher.test(contingency_table)
  
  # Store p-value and odds ratio
  p_values <- c(p_values, fisher_test$p.value)
  odds_ratios <- c(odds_ratios, fisher_test$estimate)  # Odds ratio
  
  # Store results
  fisher_results[[gene]] <- list(
    table = contingency_table,
    p_value = fisher_test$p.value,
    odds_ratio = fisher_test$estimate
  )
}

# Apply Holm Correction
adjusted_p_values <- p.adjust(p_values, method = "holm")

# Map corrected p-values back to genes
corrected_results <- data.frame(
  Gene = names(fisher_results),
  Raw_P_Value = p_values,
  Holm_Corrected_P = adjusted_p_values,
  Odds_Ratio = odds_ratios
)

# Gene Raw_P_Value Holm_Corrected_P Odds_Ratio
# 1   ABL1   0.7358377        1.0000000   1.482840
# 2  DOCK4   0.5535308        1.0000000   1.462020
# 3   XKR3   0.1215725        0.6078625   2.162200
# 4 IMMP2L   1.0000000        1.0000000   1.465570
# 5   XACT   0.2879676        1.0000000   2.237945

## apparently this wasn't the goal, need to test
## all of the relaxed filtering criteria...

# fill in the gene_name
library(dplyr)
library(stringr)

storm_k562_sc_tes.quant.filt.relaxed.clean <- storm_k562_sc_tes.quant.filt.relaxed %>%
  left_join(storm_k562_sc_tes %>% select(transcript_id, tss_gene_name), by = "transcript_id") %>%
  mutate(gene_name = if_else(tss_gene_name != ".", tss_gene_name, transcript_id)) %>%
  select(-tss_gene_name)

# aggregate stringtie_tpm per gene_name and sample
storm_k562_sc_tes_relax_aggregated_data <- storm_k562_sc_tes.quant.filt.relaxed.clean %>%
  group_by(gene_name, remapped) %>%
  summarise(sum_tpm = sum(stringtie_tpm, na.rm=TRUE), .groups="drop")  # Sum over all transcript_id

# make the expression matrix
storm_k562_sc_tes_relax_expression_matrix <- storm_k562_sc_tes_relax_aggregated_data %>%
  spread(key=remapped, value=sum_tpm, fill=0)  # Convert to wide format

# storm_k562_sc_tes_relax_expression_matrix.log <- storm_k562_sc_tes.quant.filt.relaxed.clean %>%
#   select(transcript_id, remapped, stringtie_tpm) %>%
#   mutate(log2_tpm = log2(stringtie_tpm + 1)) %>% 
#   pivot_wider(id_cols = transcript_id,  # Use only transcript_id as identifier
#               names_from = remapped, 
#               values_from = log2_tpm, 
#               values_fill = list(log2_tpm = 0))
# storm_k562_sc_tes_relax_expression_matrix.log <- as.data.frame(storm_k562_sc_tes_relax_expression_matrix.log)
# rownames(storm_k562_sc_tes_relax_expression_matrix.log) <- storm_k562_sc_tes_relax_expression_matrix.log$transcript_id
# storm_k562_sc_tes_relax_expression_matrix.log$transcript_id <- NULL

storm_k562_sc_tes_relax_expression_matrix <- as.data.frame(storm_k562_sc_tes_relax_expression_matrix)
rownames(storm_k562_sc_tes_relax_expression_matrix) <- storm_k562_sc_tes_relax_expression_matrix$gene_name
storm_k562_sc_tes_relax_expression_matrix$gene_name <- NULL
storm_k562_sc_tes_relax_expression_matrix.log <- log2(storm_k562_sc_tes_relax_expression_matrix + 1)

exp_relax.filt <- colnames(storm_k562_sc_tes_relax_expression_matrix.log) %in% binary_mat_fusions$remapped
# FALSE  TRUE 
# 9    87

# filter to common cells
storm_k562_sc_tes_relax_expression_matrix.log.filt <- storm_k562_sc_tes_relax_expression_matrix.log[,exp_relax.filt]
binary_mat_fusions.filt <- binary_mat_fusions[binary_mat_fusions$remapped %in% colnames(storm_k562_sc_tes_relax_expression_matrix.log.filt),]

all(colnames(storm_k562_sc_tes_relax_expression_matrix.log.filt) == binary_mat_fusions.filt$remapped)
# FALSE
all(colnames(storm_k562_sc_tes_relax_expression_matrix.log.filt) %in% binary_mat_fusions.filt$remapped)
# TRUE
# remap
fus.m <- match(colnames(storm_k562_sc_tes_relax_expression_matrix.log.filt),
               binary_mat_fusions$remapped)
binary_mat_fusions.filt <- binary_mat_fusions[fus.m,]
all(colnames(storm_k562_sc_tes_relax_expression_matrix.log.filt) == binary_mat_fusions.filt$remapped)
# TRUE

# recheck now that we've filtered away a bunch of cells that it still meets
# the 20-90% expression criteria
gene_nonzero_proportion <- rowSums(storm_k562_sc_tes_relax_expression_matrix.log.filt != 0) / 
  ncol(storm_k562_sc_tes_relax_expression_matrix.log.filt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1839  0.2644  0.3678  0.4332  0.5747  1.0000 
# close enough.
range_filter <- gene_nonzero_proportion >= 0.2 & gene_nonzero_proportion <= 0.8
# range_filter <- gene_nonzero_proportion >= 0.3
storm_k562_sc_tes_relax_expression_matrix.log.refilt <- storm_k562_sc_tes_relax_expression_matrix.log.filt[range_filter,]

# rip through and do wilcoxon rank sum test
library(tibble)

# expression filter?
# storm_k562_sc_tes_relax_expression_matrix.log.refilt.exp_filt <- storm_k562_sc_tes_relax_expression_matrix.log.refilt[rowMeans(as.matrix(storm_k562_sc_tes_relax_expression_matrix.log.refilt)) >= 2,]

# Convert row names to a column for easier joining
expr_df <- as.data.frame(storm_k562_sc_tes_relax_expression_matrix.log.refilt) %>%
  rownames_to_column(var = "gene_name")

# Ensure the fusion_df has sample names as a column
fusion_df <- binary_mat_fusions.filt
colnames(fusion_df)[1] <- "sample"

# Prepare an empty matrix to store p-values
fusion_genes <- colnames(fusion_df)[-1]  # Exclude the sample name column
gene_names <- expr_df$gene_name

# fisher's exact
p_values_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(fusion_genes),
                          dimnames = list(gene_names, fusion_genes))

odds_ratio_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(fusion_genes),
                            dimnames = list(gene_names, fusion_genes))

# wilcoxon rank sum
# library(dplyr)
# library(rstatix)
# library(coin)
# Prepare matrices to store p-values and effect sizes
# p_values_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(fusion_genes),
#                           dimnames = list(gene_names, fusion_genes))
# effect_size_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(fusion_genes),
#                              dimnames = list(gene_names, fusion_genes))

# Reshape expression data to long format for merging
expr_long <- expr_df %>%
  pivot_longer(cols = -gene_name, names_to = "sample", values_to = "expression")

# Loop through each fusion gene
# fisher's exact
for (fusion in fusion_genes) {

  # Get fusion status (binary: 0 or 1) for each sample
  fusion_status <- fusion_df %>% select(sample, all_of(fusion))

  # Merge expression data with fusion status
  merged_data <- expr_long %>%
    inner_join(fusion_status, by = "sample") %>%
    mutate(expression_binary = ifelse(expression > 0, 1, 0))  # Convert to expressed (1) / not expressed (0)

  # Loop through each gene
  for (gene in gene_names) {
    subset_data <- merged_data %>% filter(gene_name == gene)

    # Construct contingency table
    contingency_table <- table(subset_data$expression_binary, subset_data[[fusion]])

    # Ensure the table has at least two rows and two columns (for Fisher's test)
    if (all(dim(contingency_table) == c(2, 2))) {
      test_result <- fisher.test(contingency_table)
      p_values_matrix[gene, fusion] <- test_result$p.value
      odds_ratio_matrix[gene, fusion] <- test_result$estimate  # Extract odds ratio
    }
  }
}

# wilcoxon rank sum
# for (fusion in fusion_genes) {
# 
#   # Get fusion status (ensuring it's a factor) for each sample
#   fusion_status <- fusion_df %>% 
#     select(sample, all_of(fusion)) %>% 
#     mutate(across(all_of(fusion), as.factor))
#   
#   # Merge expression data with fusion status
#   merged_data <- expr_long %>%
#     inner_join(fusion_status, by = "sample")
#   
#   # Loop through each gene
#   for (gene in gene_names) {
#     subset_data <- merged_data %>% filter(gene_name == gene)
#     
#     # Skip if the fusion status doesn't vary
#     if (length(unique(dplyr::pull(subset_data, fusion))) < 2) next
#     
#     # Mutate the fusion column to a standard name, e.g. "fusion_status"
#     subset_data <- subset_data %>% mutate(fusion_status = dplyr::pull(subset_data, fusion))
#     
#     # Run the Wilcoxon rank sum test on the continuous expression values
#     test_result <- wilcox.test(expression ~ fusion_status, data = subset_data)
#     p_values_matrix[gene, fusion] <- test_result$p.value
#     
#     # Compute the effect size (e.g., rank-biserial correlation) using rstatix
#     eff_size <- subset_data %>% wilcox_effsize(expression ~ fusion_status)
#     effect_size_matrix[gene, fusion] <- eff_size$effsize
#   }
# }

# logistic
# for (fusion in fusion_genes) {
#   
#   # Get fusion status (binary: 0 or 1) for each sample
#   fusion_status <- fusion_df %>% 
#     select(sample, all_of(fusion)) %>% 
#     mutate(across(all_of(fusion), factor))
# 
#   # Merge expression data with fusion status
#   merged_data <- expr_long %>%
#     inner_join(fusion_status, by = "sample")
#   
#   # Loop through each gene
#   for (gene in gene_names) {
#     subset_data <- merged_data %>% filter(gene_name == gene)
#     
#     # Only run logistic regression if the fusion outcome varies
#     if (length(unique(subset_data[[fusion]])) > 1) {
#       # Fit logistic regression: outcome is fusion status and predictor is gene expression
#       fusion_var <- paste0("`", fusion, "`")
#       form <- as.formula(paste(fusion_var, "~ expression"))
#       model <- glm(form, data = subset_data, family = "binomial")
#       # Extract the coefficient summary
#       coefs <- summary(model)$coefficients
#       
#       # Save the p-value and the odds ratio (exponentiated coefficient for 'expression')
#       p_values_matrix[gene, fusion] <- coefs["expression", "Pr(>|z|)"]
#       odds_ratio_matrix[gene, fusion] <- exp(coefs["expression", "Estimate"])
#     }
#   }
# }

# Convert matrices to data frames
p_values_df <- as.data.frame(p_values_matrix)
odds_ratios_df <- as.data.frame(odds_ratio_matrix)
# effect_size_df <- as.data.frame(effect_size_matrix)

pvals_and_oddratios <- cbind(p_values_df[,c(1:2)],
                             odds_ratios_df[,c(1:2)])
colnames(pvals_and_oddratios) <- c("BCR_ABL1.pval",
                                   "NUP214_XKR3.pval",
                                   "BCR_ABL1.or",
                                   "NUP214_XKR3.or")
# pvals_and_oddratios$bcrabl_padj <- p.adjust(pvals_and_oddratios$BCR_ABL1.pval,
#                                             method = "BH")

# pvals_and_effsizes <- cbind(p_values_df[,c(1:2)],
#                             effect_size_df[,c(1:2)])
# colnames(pvals_and_effsizes) <- c("BCR_ABL1.pval",
#                                    "NUP214_XKR3.pval",
#                                    "BCR_ABL1.effsize",
#                                    "NUP214_XKR3.effsize")
# pvals_and_effsizes$bcrabl_padj <- p.adjust(pvals_and_effsizes$BCR_ABL1.pval,
#                                            method = "BH")
# pvals_and_effsizes$nupxkr_padj <- p.adjust(pvals_and_effsizes$NUP214_XKR3.pval,
#                                            method = "BH")
# BCR-ABL1
bcr_abl1_te_candidates <- pvals_and_oddratios[,c(1,3)]
bcr_abl1_te_candidates <- bcr_abl1_te_candidates[bcr_abl1_te_candidates$BCR_ABL1.pval < 0.05,]
bcr_abl1_te_candidates <- bcr_abl1_te_candidates[bcr_abl1_te_candidates$BCR_ABL1.or >= 3,]

# TPM filter
# bcr_abl1_te_candidates <- bcr_abl1_te_candidates[rowSums(bcr_abl1_te_candidates) >= 3,]

# back map the genes
# bcr_txps <- rownames(bcr_abl1_te_candidates)
# bcr_txps.sc_tes <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% bcr_txps,]
# bcr_txps.sc_tes.tu <- bcr_txps.sc_tes$transcript_id
# bcr_txps.sc_tes.annot <- ifelse(bcr_txps.sc_tes$tss_gene_name != ".",
#                                 bcr_txps.sc_tes$tss_gene_name,
#                                 bcr_txps.sc_tes$transcript_id)
# names(bcr_txps.sc_tes.tu) <- bcr_txps.sc_tes.annot

# add the genes (if known) and convert to syntax of AluJ-TP53 for row names
bcr_abl1_te_candidates.genes.gtf <- gtf.genes[gtf.genes$gene_id %in% rownames(bcr_abl1_te_candidates),]
gtf.genes.bcr_abl_ens_ids <- bcr_abl1_te_candidates.genes.gtf$gene_name
names(gtf.genes.bcr_abl_ens_ids) <- bcr_abl1_te_candidates.genes.gtf$gene_id

bcr_abl1_te_candidates.genes.ens <- ifelse(rownames(bcr_abl1_te_candidates) %in% names(gtf.genes.bcr_abl_ens_ids),
                                       gtf.genes.bcr_abl_ens_ids[rownames(bcr_abl1_te_candidates)],
                                       rownames(bcr_abl1_te_candidates))
names(bcr_abl1_te_candidates.genes.ens) <- rownames(bcr_abl1_te_candidates)

# look up to get the TE family part
# ensembl genes
storm_k562_sc_tes.bcr_abl_candidates.ens_id <- storm_k562_sc_tes[storm_k562_sc_tes$tss_gene_name %in% names(bcr_abl1_te_candidates.genes.ens),]
storm_k562_sc_tes.bcr_abl_candidates.tu_id <- storm_k562_sc_tes.bcr_abl_candidates.ens_id[storm_k562_sc_tes.bcr_abl_candidates.ens_id$transcript_id %in% storm_k562_sc_tes.quant.filt.relaxed.clean$transcript_id,]

# resolve duplicate genes with primary TE-txp
# ENSG00000167562 is primarily driven by TU6547/AluYc
# ENSG00000203727 is primarily driven by TU13092/MSTA
# ENSG00000257354 is primarily driven by TU2947/MLT1G
# ENSG00000286403 is primarily driven by TU11938/LTR12

bcr_abl_ens_te <- storm_k562_sc_tes.bcr_abl_candidates.tu_id$TE_subfamily
names(bcr_abl_ens_te) <- storm_k562_sc_tes.bcr_abl_candidates.tu_id$tss_gene_name

primary_te <- c(
  "ENSG00000167562" = "AluYc",
  "ENSG00000203727" = "MSTA",
  "ENSG00000257354" = "MLT1G",
  "ENSG00000286403" = "LTR12"
)

# convert your named vector to a df/tibble
bcr_abl_ens_te.df <- tibble(
  gene = names(bcr_abl_ens_te),
  te = bcr_abl_ens_te
)

# add the primary TE (if available) for each gene
bcr_abl_ens_te.df <- bcr_abl_ens_te.df %>%
  mutate(primary = coalesce(primary_te[gene], te))

bcr_abl_ens_te.df_resolved <- bcr_abl_ens_te.df %>%
  # Keep only rows where the original TE equals the primary TE
  filter(te == primary) %>%
  # If there are duplicate rows with the same gene, te, and primary, keep only the first
  distinct(gene, te, primary, .keep_all = TRUE)

# convert back to a named vector
bcr_abl_ens_te <- bcr_abl_ens_te.df_resolved$te
names(bcr_abl_ens_te) <- bcr_abl_ens_te.df_resolved$gene

# ENSG00000126091 ENSG00000143106 ENSG00000266028 ENSG00000110696 ENSG00000149716 
# "MIRb"          "HAL1"         "AluJr"         "AluSp"         "MLT1F" 
# ENSG00000256537 ENSG00000257354 ENSG00000186908 ENSG00000108826 ENSG00000121060 
# "MIRb"         "MLT1G"         "HAL1b"        "AluSq2"           "MIR" 
# ENSG00000167562 ENSG00000286403 ENSG00000203727 ENSG00000146707 ENSG00000174469 
# "AluYc"         "LTR12"          "MSTA"        "MER41A"         "L1PA5" 
# ENSG00000130119 
# "AluSg" 

# TEProf3 TE transcript IDs
storm_k562_sc_tes.bcr_abl_candidates.tu <- storm_k562_sc_tes[storm_k562_sc_tes$transcript_id %in% rownames(bcr_abl1_te_candidates),]
bcr_abl_tu_te <- storm_k562_sc_tes.bcr_abl_candidates.tu$TE_subfamily
names(bcr_abl_tu_te) <- storm_k562_sc_tes.bcr_abl_candidates.tu$transcript_id
# TU627      TU1245      TU2602      TU2754      TU2966      TU3123 
# "AluY" "HERVH-int"     "L1MD1"     "L1MEb"     "L1MB7"       "L2a" 
# TU5616      TU9146      TU9136     TU12551     TU14281     TU14253 
# "AluY"    "MER61A"    "AluYb8"     "AluJb"    "LTR78B"      "L1M1" 
# TU14717 
# "LTR12C" 

bcr_abl_te_subfam <- c(bcr_abl_ens_te,
                       bcr_abl_tu_te)

bcr_abl_te_subfam.clean <- sapply(names(bcr_abl_te_subfam), function(gene_id) {
  # Look up the gene name
  gene_name <- gtf.genes.bcr_abl_ens_ids[gene_id]
  if (!is.na(gene_name)) {
    # For ENSG entries: use TE-GeneName format
    paste(bcr_abl_te_subfam[gene_id], gene_name, sep = "-")
  } else {
    # For TU (or any gene_id not found): use GeneID-TE format
    paste(bcr_abl_te_subfam[gene_id], gene_id, sep = "-")
  }
})

# NUP214-XKR3
nup214_xkr3_te_candidates <- pvals_and_oddratios[,c(2,4)]
nup214_xkr3_te_candidates <- nup214_xkr3_te_candidates[nup214_xkr3_te_candidates$NUP214_XKR3.pval < 0.05,]
nup214_xkr3_te_candidates <- nup214_xkr3_te_candidates[nup214_xkr3_te_candidates$NUP214_XKR3.or >= 3,]

# TPM filter
nup214_xkr3_te_candidates <- nup214_xkr3_te_candidates[rowSums(nup214_xkr3_te_candidates) >= 3,]

# any common TE-txps or genes?
table(rownames(bcr_abl1_te_candidates) %in% rownames(nup214_xkr3_te_candidates))
# FALSE  TRUE 
# 27     1 

# what are they?
common_te_txps <- rownames(bcr_abl1_te_candidates[rownames(bcr_abl1_te_candidates) %in% rownames(nup214_xkr3_te_candidates),])
# [1] "ENSG00000143106"
# proteosome subunit - PSMA5

# let's put a heatmap together with the BCR-ABL1 positive cells as annotation
bcr_abl_annots <- factor(binary_mat_fusions.filt$`BCR--ABL1`)
names(bcr_abl_annots) <- ifelse(bcr_abl_annots == 0, "YES", "NO")
bcr_abl_te_txps <- storm_k562_sc_tes_relax_expression_matrix.log.filt[rownames(storm_k562_sc_tes_relax_expression_matrix.log.filt) %in% 
                                                                        rownames(bcr_abl1_te_candidates),]
bcr_abl_te_txps.rnames <- rownames(bcr_abl_te_txps)
# Replace those that match names in combined_vector with the corresponding new value
bcr_abl_te_txps.idx <- bcr_abl_te_txps.rnames %in% names(bcr_abl_te_subfam.clean)
bcr_abl_te_txps.rnames[bcr_abl_te_txps.idx] <- bcr_abl_te_subfam.clean[bcr_abl_te_txps.rnames[bcr_abl_te_txps.idx]]
rownames(bcr_abl_te_txps) <- bcr_abl_te_txps.rnames

# let's put a heatmap together with the NUP214-XKR3 positive cells as annotation
nup214_xkr3_annots <- factor(binary_mat_fusions.filt$`NUP214--XKR3`)
names(nup214_xkr3_annots) <- ifelse(nup214_xkr3_annots == 0, "YES", "NO")
nup214_xkr3_te_txps <- storm_k562_sc_tes_relax_expression_matrix.log.filt[rownames(storm_k562_sc_tes_relax_expression_matrix.log.filt) %in% 
                                                                            rownames(nup214_xkr3_te_candidates),]


# helper
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

# expression helper func
getColorRamp <- function(low = "#FFFFFF", high = "darkred") {
  return(circlize::colorRamp2(seq(0, 6, length.out = 2), c(low, high)))
}

# filter to only common cells for plotting
# storm_sce.filt.sub <- storm_sce.filt[,colnames(storm_sce.filt) %in% gsub("_teprof3", "", colnames(bcr_abl_te_txps))]

# re-org
# cell.m <- match(gsub("_teprof3", "", colnames(bcr_abl_te_txps)),
#                 colnames(storm_sce.filt.sub))
# storm_sce.filt.sub <- storm_sce.filt.sub[,cell.m]

#careful
# all(colnames(storm_sce.filt.sub) == gsub("_teprof3", "", colnames(bcr_abl_te_txps)))
# # TRUE
# 
# all(colnames(storm_sce.filt.sub) == gsub("_teprof3", "", colnames(nup214_xkr3_te_txps)))
# # TRUE

# BCR-ABL1
library(ComplexHeatmap)
hist.anno <- HeatmapAnnotation(`BCR--ABL1` = bcr_abl_annots,
                               col = list(`BCR--ABL1` = c("1" = "black",
                                                          "0" = "gray")),
                               annotation_legend_param = list(
                                 `BCR--ABL1` = list(
                                   title = "BCR--ABL1",
                                   at = levels(bcr_abl_annots),
                                   labels = c("NO", "YES"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(4, "cm"))))
# filter out AluSp-C11orf58
bcr_abl_te_txps.filt <- bcr_abl_te_txps[!rownames(bcr_abl_te_txps) %in% "AluSp-C11orf58",]
annot_filt <- !colnames(bcr_abl_te_txps.filt) %in% c("K17_teprof3",
                                                     "F3_teprof3")
bcr_abl_te_txps.filt <- bcr_abl_te_txps.filt[,!colnames(bcr_abl_te_txps.filt) %in% c("K17_teprof3",
                                                                                     "F3_teprof3")]
bcr_abl_annots <- bcr_abl_annots[annot_filt]

heatmap_legend_params <- list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                              labels_gp = gpar(fontsize = 12),
                              legend_width = unit(4, "cm"),
                              legend_height = unit(4, "cm"))

Heatmap(as.matrix(bcr_abl_te_txps.filt),
        top_annotation = hist.anno,
        column_split = bcr_abl_annots,
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_title = NULL,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

# great - now overlay the MYC and PCNA expression profiles as additional annotations
library(scran)
storm_sce.filt <- computeSpikeFactors(storm_sce.filt, "ERCC")
storm_sce.filt <- computeSumFactors(storm_sce.filt)
storm_sce.filt.norm <- logNormCounts(storm_sce.filt)

# MYC
myc_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000136997",]

# PCNA
pcna_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000132646",]

# NPM1 as a proxy for ribosome assembly
npm1_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000107833",]

# RUNX2
runx2_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000124813",]

# MECOM
mecom_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000085276",]

# CD20
cd20_storm <- logcounts(storm_sce.filt.norm)[rownames(storm_sce.filt.norm) %in% "ENSG00000156738",]
# NOT EXPRESSED


# filter to the cells in the te exp matrix
sc_filt <- paste0(colnames(storm_sce.filt.norm),
                  "_teprof3") %in% colnames(bcr_abl_te_txps.filt)
myc_storm <- myc_storm[sc_filt]
pcna_storm <- pcna_storm[sc_filt]
npm1_storm <- npm1_storm[sc_filt]

myc_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                     paste0(names(myc_storm),
                            "_teprof3"))
myc_storm <- myc_storm[myc_storm.m]
names(myc_storm) <- paste0(names(myc_storm),
                           "_teprof3")

pcna_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                     paste0(names(pcna_storm),
                            "_teprof3"))
pcna_storm <- pcna_storm[pcna_storm.m]
names(pcna_storm) <- paste0(names(pcna_storm),
                           "_teprof3")

npm1_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                      paste0(names(npm1_storm),
                             "_teprof3"))
npm1_storm <- npm1_storm[npm1_storm.m]
names(npm1_storm) <- paste0(names(npm1_storm),
                            "_teprof3")

mecom_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                      paste0(names(mecom_storm),
                             "_teprof3"))
mecom_storm <- mecom_storm[mecom_storm.m]
names(mecom_storm) <- paste0(names(mecom_storm),
                            "_teprof3")

runx2_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                       paste0(names(runx2_storm),
                              "_teprof3"))
runx2_storm <- runx2_storm[runx2_storm.m]
names(runx2_storm) <- paste0(names(runx2_storm),
                             "_teprof3")

cd20_storm.m <- match(colnames(bcr_abl_te_txps.filt),
                       paste0(names(cd20_storm),
                              "_teprof3"))
cd20_storm <- cd20_storm[cd20_storm.m]
names(cd20_storm) <- paste0(names(cd20_storm),
                             "_teprof3")

all(names(pcna_storm) == colnames(bcr_abl_te_txps.filt))
all(names(myc_storm) == colnames(bcr_abl_te_txps.filt))

# do something dumb and convert it to % of total expression
pcna_storm.scale <- scales::rescale(pcna_storm, c(0,1))
myc_storm.scale <- scales::rescale(myc_storm, c(0,1))

# try ordering based on myc expression
myc_storm.ord <- myc_storm[order(myc_storm, decreasing = T)]
pcna_storm.ord.m <- match(names(myc_storm.ord),
                          names(pcna_storm))
pcna_storm.ord <- pcna_storm[pcna_storm.ord.m]
bcr_abl_annots <- bcr_abl_annots[pcna_storm.ord.m]
bcr_abl_te_txps.filt <- bcr_abl_te_txps.filt[,pcna_storm.ord.m]

# replot
hist.anno <- HeatmapAnnotation(`BCR--ABL1` = bcr_abl_annots,
                               myc_expression = myc_storm,
                               pcna_expression = pcna_storm,
                               npm3_expression = npm1_storm,
                               mecom_expression = mecom_storm,
                               runx2_expression = runx2_storm,
                               col = list(`BCR--ABL1` = c("1" = "black",
                                                          "0" = "gray"),
                                          myc_expression = getColorRamp(),
                                          pcna_expression = getColorRamp(),
                                          npm3_expression = getColorRamp(),
                                          mecom_expression = getColorRamp(),
                                          runx2_expression = getColorRamp()),
                               annotation_legend_param = list(
                                 `BCR--ABL1` = list(
                                   title = "BCR--ABL1",
                                   at = levels(bcr_abl_annots),
                                   labels = c("NO", "YES"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(4, "cm")),
                                 myc_expression = list(
                                   title = "MYC",
                                   at = c(0, 6),
                                   labels = c("0", "6"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(2, "cm")),
                                 pcna_expression = list(
                                   title = "PCNA",
                                   at = c(0, 6),
                                   labels = c("0", "6"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(2, "cm")),
                                 npm3_expression = list(
                                   title = "NPM3",
                                   at = c(0, 6),
                                   labels = c("0", "6"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(2, "cm")),
                                 mecom_expression = list(
                                   title = "MECOM",
                                   at = c(0, 6),
                                   labels = c("0", "6"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(2, "cm")),
                                 runx2_expression = list(
                                   title = "RUNX2",
                                   at = c(0, 6),
                                   labels = c("0", "6"),
                                   title_gp = gpar(fontsize = 12, fontface = "bold"),
                                   labels_gp = gpar(fontsize = 12),
                                   legend_width = unit(4, "cm"),
                                   legend_height = unit(2, "cm"))))
Heatmap(as.matrix(bcr_abl_te_txps.filt),
        top_annotation = hist.anno,
        column_split = bcr_abl_annots,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        column_title = NULL,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)",
        heatmap_legend_param = heatmap_legend_params)

library(ggplot2)
to_plot <- data.frame(bcrabl = bcr_abl_annots,
                      myc = myc_storm,
                      pcna = pcna_storm)

ggplot(to_plot, aes(x = myc, y = pcna, color = bcr_abl_annots)) +
  geom_point() +
  theme_bw(12)

# NUP214-XKR3
library(ComplexHeatmap)
hist.anno <- HeatmapAnnotation(`NUP214--XKR3` = nup214_xkr3_annots,
                               col = list(`NUP214--XKR3` = c("1" = "black",
                                                          "0" = "gray")),
                               annotation_legend_param = list(
                                 `NUP214--XKR3` = list(
                                   title = "NUP214--XKR3",
                                   at = levels(nup214_xkr3_annots),
                                   labels = c("NO", "YES"))))

Heatmap(as.matrix(nup214_xkr3_te_txps),
        top_annotation = hist.anno,
        column_split = nup214_xkr3_annots,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = NULL,
        col = getJetColors(circlize.cols = TRUE),
        name = "log2(TPM+1)")



# Gviz view of the 




