# restructure the bed files
# to convert to a GTF

library(rtracklayer)

# unidirectional known enhancers
uni_enh <- read.delim("crispr_mpra_starrseq_grocap_unidirectional_1kb_filt_1kb_resize.hg38.distal_only.bed.gz",
                      header = FALSE)
colnames(uni_enh) <- c("chr",
                       "start",
                       "end",
                       "name",
                       "score",
                       "strand")
# divergent known enhancers
div_enh <- read.delim("crispr_mpra_starrseq_grocap_divergent_1kb_filt_1kb_resize.hg38.distal_only.bed.gz",
                      header = FALSE)
colnames(div_enh) <- c("chr",
                       "start",
                       "end",
                       "name",
                       "score",
                       "strand")

# reformat to GTF
uni_enh.gr <- as(uni_enh, "GRanges")
div_enh.gr <- as(div_enh, "GRanges")

# assign specific names so we know which is which
uni_enh.gr$name <- paste0("enh_uni_", uni_enh.gr$name)
div_enh.gr$name <- paste0("enh_div_", div_enh.gr$name)

# combine
all_enh.gr <- c(uni_enh.gr,
                div_enh.gr)
all_enh.gr <- sort(all_enh.gr)

all_enh.df <- as.data.frame(all_enh.gr)

# convert and annotate
create_gtf <- function(df) {
  # Create gene rows
  gene_rows <- data.frame(
    seqname = df$seqnames,
    source = "PINTS",
    feature = "gene",
    start = df$start,
    end = df$end,
    score = 0.000000,
    strand = df$strand,
    frame = ".",
    attributes = paste0(
      'gene_id "', df$name, '"; transcript_id "', df$name,
      '"; gene_type "enhancer_annot"; gene_status "KNOWN"; ',
      'gene_name "', df$name, '"; transcript_type "enhancer_annot"; ',
      'transcript_status "KNOWN"; transcript_name "', df$name, '"; level 2;'
    )
  )
  
  # Create transcript rows
  transcript_rows <- data.frame(
    seqname = df$seqnames,
    source = "PINTS",
    feature = "transcript",
    start = df$start,
    end = df$end,
    score = 0.000000,
    strand = df$strand,
    frame = ".",
    attributes = paste0(
      'gene_id "', df$name, '"; transcript_id "', df$name,
      '"; gene_type "enhancer_annot"; gene_status "KNOWN"; ',
      'gene_name "', df$name, '"; transcript_type "enhancer_annot"; ',
      'transcript_status "KNOWN"; transcript_name "', df$name, '"; level 2;'
    )
  )
  
  # Create exon rows
  exon_rows <- data.frame(
    seqname = df$seqnames,
    source = "PINTS",
    feature = "exon",
    start = df$start,
    end = df$end,
    score = 0.000000,
    strand = df$strand,
    frame = ".",
    attributes = paste0(
      'gene_id "', df$name, '"; transcript_id "', df$name, '";'
    )
  )
  
  # Combine rows in the desired order: gene, transcript, exon
  gtf <- rbind(gene_rows, transcript_rows, exon_rows)
  
  return(gtf)
}

all_enh.df.gtf <- create_gtf(all_enh.df)

# Write to a GTF file
write.table(
  all_enh.df.gtf,
  file = "crispr_mpra_starrseq_grocap_uni_div_1kb_filt_1kb_resize.hg38.distal_only.gtf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# export
# rtracklayer::export(all_enh.gr, "crispr_mpra_starrseq_grocap_uni_div_1kb_filt_1kb_resize.hg38.distal_only.gtf")

# rebuild for the pints distal enhancer set for K562 as well
pints_distal_k562 <- read.delim("~/Documents/manuscripts/storm_seq/erna/Distal_K562_gro_pro_cap_pints_hg38_gencode24.bed",
                                header = FALSE)
pints_distal_k562 <- pints_distal_k562[,c(1:3)]
colnames(pints_distal_k562) <- c("chr",
                                 "start",
                                 "end")
pints_distal_k562$name <- paste0("pints_k562_distal_enh_",
                                 seq(1, nrow(pints_distal_k562)))
pints_distal_k562$strand <- "*"
pints_distal_k562.gr <- as(pints_distal_k562, "GRanges")
seqlevelsStyle(pints_distal_k562.gr) <- "Ensembl"

# do the filtering and resizing
summary(width(pints_distal_k562.gr))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 101.0   401.0   401.0   494.4   541.0  3584.0 

# filter things >1kb
pints_distal_k562.gr.filt <- pints_distal_k562.gr[!width(pints_distal_k562.gr) > 1000,]
# resize from center
pints_distal_k562.gr.filt.resize <- resize(pints_distal_k562.gr.filt,
                                           width = 1000,
                                           fix = "center")
pints_distal_k562.gr.filt.resize <- sort(pints_distal_k562.gr.filt.resize)
pints_distal_k562.df <- as.data.frame(pints_distal_k562.gr.filt.resize)

pints_enh.df.gtf <- create_gtf(pints_distal_k562.df)

# Write to a GTF file
write.table(
  pints_enh.df.gtf,
  file = "~/Documents/manuscripts/storm_seq/erna/pints_distal_grocap_procap_1kb_filt_1kb_resize.hg38.distal_only.gtf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# write out as a bed too for profile plots
rtracklayer::export(pints_distal_k562.gr.filt.resize,
                    "~/Documents/manuscripts/storm_seq/erna/pints_distal_grocap_procap_1kb_filt_1kb_resize.hg38.distal_only.bed")

## ingest the starsolo results
library(Matrix)
library(SingleCellExperiment)

# write a small function to annotate the matrix market format
read_starsolo_gene <- function(res_dir, cell_name=NULL) {
  ## list the files and find the mtx
  mtx <- list.files(res_dir,
                    pattern = "UniqueAndMult-EM.mtx",
                    full.names = TRUE)
  # read in
  mat.mm <- readMM(mtx)
  mat.mm <- as(mat.mm, "dgCMatrix")
  genes <- list.files(res_dir,
                      pattern = "features.tsv",
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

storm_counts_files <- list.files("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/full_depth/k562",
                                 pattern = "UniqueAndMult-EM.mtx",
                                 recursive = TRUE,
                                 full.names = TRUE)
storm_counts_files <- dirname(storm_counts_files)

first_gsub <- gsub("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/full_depth/k562/",
                   "", storm_counts_files)
names(storm_counts_files) <- gsub("_solooutput/GeneFull/raw",
                                  "", first_gsub)

# read in data
storm_counts_raw <- lapply(names(storm_counts_files), function(x) {
  # get the path
  kb_mm <- read_starsolo_gene(storm_counts_files[x],
                              x)
  return(kb_mm)
})

# merge
storm_counts_raw.mat <- do.call(cbind, storm_counts_raw)

# make the SCE object
storm_sce.raw <- SingleCellExperiment(assays = list(counts = storm_counts_raw.mat))

# get the detected genes
# should round these...
storm_sce.raw$NumGenesExpressed <- colSums2(counts(storm_sce.raw) > 0)

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

# fix
storm_sce.raw <- fix_rownames(storm_sce.raw)
# Fixed 60671 rownames.

# pull out the ERCCs
ercc.sce <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "ERCC") <- ercc.sce
storm_sce.raw <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw), invert = TRUE),]

# create a separate altexp for the enhancers
enh.sce <- storm_sce.raw[grep("enh_", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "Enhancers") <- enh.sce
storm_sce.raw <- storm_sce.raw[grep("enh_", rownames(storm_sce.raw), invert = TRUE),]

# storm_sce.raw <- readRDS("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_raw_sce_with_known_k562_enh_annots_hg38_ens101.rds")

# QC
library(scuttle)
storm_sce.raw <- addPerCellQCMetrics(storm_sce.raw)
filt_reasons <- perCellQCFilters(storm_sce.raw,
                                 sub.fields = c("altexps_ERCC_percent"))
colSums(as.matrix(filt_reasons))
# low_lib_size            low_n_features high_altexps_ERCC_percent 
# 2                         4                         7 
# discard 
# 9

# filter
storm_sce.filt <- storm_sce.raw[,!filt_reasons$discard]

saveRDS(storm_sce.filt,
        file = "/Users/ben.johnson/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_filt_sce_with_known_k562_enh_annots_hg38_ens101.rds")
#storm_sce.filt <- readRDS("/Volumes/projects/shen/projects/2021_02_25_ben_projects/laptop_backup/kb_python/kb_python_storm_test/storm_100k_sce_filt_hg38_ens101.rds")

# how many known enhancers per cell
storm_sce.filt$NumEnhExpressed <- colSums2(counts(altExp(storm_sce.filt, "Enhancers")) > 0)
# re-calculate num genes
storm_sce.filt$NumGenesExpressed <- colSums2(counts(storm_sce.filt) > 0)
# re-calculate depth
storm_sce.filt$depth_size <- colSums2(counts(storm_sce.filt))

library(ggplot2)
library(viridis)

ggplot(colData(storm_sce.filt), aes(x = NumEnhExpressed, y = depth_size, color = depth_size)) +
  geom_point() +
  scale_color_viridis_c(end = 0.8, option = "magma") +
  theme_bw(12)

## vasa
vasa_counts_files <- list.files("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/vasa_plate/full_depth",
                                 pattern = "UniqueAndMult-EM.mtx",
                                 recursive = TRUE,
                                 full.names = TRUE)
vasa_counts_files <- dirname(vasa_counts_files)

first_gsub <- gsub("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/vasa_plate/full_depth/",
                   "", vasa_counts_files)
names(vasa_counts_files) <- gsub("output/GeneFull/raw",
                                  "", first_gsub)

# read in data
vasa_counts_raw <- lapply(names(vasa_counts_files), function(x) {
  # get the path
  kb_mm <- read_starsolo_gene(vasa_counts_files[x],
                              x)
  return(kb_mm)
})

# merge
vasa_counts_raw.mat <- do.call(cbind, vasa_counts_raw)

# make the SCE object
vasa_sce.raw <- SingleCellExperiment(assays = list(counts = vasa_counts_raw.mat))

# get the detected genes
vasa_sce.raw$NumGenesExpressed <- colSums2(counts(vasa_sce.raw) > 0)

# pull out the ERCCs
ercc.sce <- vasa_sce.raw[grep("ERCC-", rownames(vasa_sce.raw)),]
altExp(vasa_sce.raw, "ERCC") <- ercc.sce
vasa_sce.raw <- vasa_sce.raw[grep("ERCC-", rownames(vasa_sce.raw), invert = TRUE),]

# create a separate altexp for the enhancers
enh.sce <- vasa_sce.raw[grep("enh_", rownames(vasa_sce.raw)),]
altExp(vasa_sce.raw, "Enhancers") <- enh.sce
vasa_sce.raw <- vasa_sce.raw[grep("enh_", rownames(vasa_sce.raw), invert = TRUE),]

# vasa_sce.raw <- readRDS("~/Documents/manuscripts/storm_seq/erna/k562_vasa/vasa_plate_raw_sce_with_known_k562_enh_annots_hg38_ens101.rds")

# QC
vasa_sce.raw <- addPerCellQCMetrics(vasa_sce.raw)
filt_reasons <- perCellQCFilters(vasa_sce.raw)
colSums(as.matrix(filt_reasons))
# low_lib_size low_n_features        discard 
# 54             68             68 

# filter
vasa_sce.filt <- vasa_sce.raw[,!filt_reasons$discard]

# remove controls too
p1_neg_controls <- c("vai_vasa_plate1_merged_009", "vai_vasa_plate1_merged_047",
                     "vai_vasa_plate1_merged_108", "vai_vasa_plate1_merged_123",
                     "vai_vasa_plate1_merged_207", "vai_vasa_plate1_merged_318",
                     "vai_vasa_plate1_merged_334", "vai_vasa_plate1_merged_379")
p2_neg_controls <- c("vai_vasa_plate2_merged_009", "vai_vasa_plate2_merged_047",
                     "vai_vasa_plate2_merged_108", "vai_vasa_plate2_merged_123",
                     "vai_vasa_plate2_merged_207", "vai_vasa_plate2_merged_318",
                     "vai_vasa_plate2_merged_334", "vai_vasa_plate2_merged_379")
vasa_sce.filt <- vasa_sce.filt[,!colnames(vasa_sce.filt) %in% c(p1_neg_controls,
                                                                p2_neg_controls)]

saveRDS(vasa_sce.filt,
        file = "/Users/ben.johnson/Documents/manuscripts/storm_seq/erna/k562_vasa/vasa_plate_filt_sce_with_known_k562_enh_annots_hg38_ens101.rds")
#vasa_sce.filt <- readRDS("/Volumes/projects/shen/projects/2021_02_25_ben_projects/laptop_backup/kb_python/kb_python_vasa_test/vasa_100k_sce_filt_hg38_ens101.rds")

# how many known enhancers per cell
vasa_sce.filt$NumEnhExpressed <- colSums2(counts(altExp(vasa_sce.filt, "Enhancers")) > 0)
# re-calculate num genes
vasa_sce.filt$NumGenesExpressed <- colSums2(counts(vasa_sce.filt) > 0)
# re-calculate depth
vasa_sce.filt$depth_size <- colSums2(counts(vasa_sce.filt))

library(ggplot2)
library(viridis)

ggplot(colData(vasa_sce.filt), aes(x = NumEnhExpressed, y = depth_size, color = depth_size)) +
  geom_point() +
  scale_color_viridis_c(end = 0.8, option = "magma") +
  theme_bw(12)

summary(vasa_sce.filt$NumEnhExpressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.000   2.000   2.699   4.000  10.000 
summary(storm_sce.filt$NumEnhExpressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   3.500   5.000   5.709   7.500  14.000

# bin and summarize proportion of cells
storm_enh <- data.frame(detected = storm_sce.filt$NumEnhExpressed,
                        technology = "STORM-seq")
vasa_enh <- data.frame(detected = vasa_sce.filt$NumEnhExpressed,
                       technology = "VASA-seq")
all_enh <- rbind(storm_enh,
                 vasa_enh)

library(dplyr)

df <- all_enh %>%
  mutate(
    threshold = case_when(
      detected == 0 ~ "0",       # Keep "0" temporarily for adjustment
      detected == 1 ~ "1",
      detected <= 5 ~ "2-5",
      detected <= 10 ~ "6-10",
      detected > 10 ~ ">10"
    )
  )

# Calculate counts for threshold 0 and others
df_counts <- df %>%
  group_by(technology, threshold) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(threshold = factor(threshold, levels = c("0", "1", "2-5", "6-10", ">10"))) %>%
  arrange(technology, threshold)

# Separate the counts for threshold 0
df_zero <- df_counts %>%
  filter(threshold == "0") %>%
  rename(count_zero = count) %>%
  select(-threshold)

# Adjust counts by excluding threshold 0
df_adjusted <- df_counts %>%
  filter(threshold != "0") %>%
  left_join(df_zero, by = "technology") %>%
  group_by(technology) %>%  # Group by technology for calculations
  mutate(
    adjusted_total = sum(count) + count_zero,          # Total includes threshold 0
    proportion = count / adjusted_total,              # Adjusted proportions
    cumulative_decrease = rev(cumsum(rev(proportion))) # Reverse cumulative proportions
  ) %>%
  ungroup()

df_adjusted <- rbind(df_adjusted,
                    data.frame(technology = "VASA-seq",
                               threshold = as.factor(">10"),
                               count = 0,
                               count_zero = 23,
                               adjusted_total = 417,
                               proportion = 0.000,
                               cumulative_decrease = 0.000))

ggplot(df_adjusted, aes(x = threshold, y = cumulative_decrease,
                        color = technology, group = technology)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0,
             linewidth = 1,
             linetype = "dashed",
             color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF")) +
  xlab("Experimentally Validated K-562\nIntergenic Enhancers (per-cell)") +
  ylab("Proportion of Total Cells") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# do it again with all PINTS gro/pro-cap data for K562, distal
storm_counts_files <- list.files("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/full_depth/k562",
                                 pattern = "UniqueAndMult-EM.mtx",
                                 recursive = TRUE,
                                 full.names = TRUE)
storm_counts_files <- dirname(storm_counts_files)

# grab just the pints dirs
storm_counts_files <- storm_counts_files[grep("pintsoutput", storm_counts_files)]

first_gsub <- gsub("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/full_depth/k562/",
                   "", storm_counts_files)
names(storm_counts_files) <- gsub("_solo_pintsoutput/GeneFull/raw",
                                  "", first_gsub)

# read in data
storm_counts_raw <- lapply(names(storm_counts_files), function(x) {
  # get the path
  kb_mm <- read_starsolo_gene(storm_counts_files[x],
                              x)
  return(kb_mm)
})

# merge
storm_counts_raw.mat <- do.call(cbind, storm_counts_raw)

# make the SCE object
storm_sce.raw <- SingleCellExperiment(assays = list(counts = storm_counts_raw.mat))

# get the detected genes
# should round these...
storm_sce.raw$NumGenesExpressed <- colSums2(counts(storm_sce.raw) > 0)

# pull out the ERCCs
ercc.sce <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "ERCC") <- ercc.sce
storm_sce.raw <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw), invert = TRUE),]

# create a separate altexp for the enhancers
enh.sce <- storm_sce.raw[grep("pints_k562_distal_enh", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "Enhancers") <- enh.sce
storm_sce.raw <- storm_sce.raw[grep("pints_k562_distal_enh", rownames(storm_sce.raw), invert = TRUE),]

# storm_sce.raw <- readRDS("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_raw_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")

# QC
library(scuttle)
storm_sce.raw <- addPerCellQCMetrics(storm_sce.raw)
filt_reasons <- perCellQCFilters(storm_sce.raw,
                                 sub.fields = c("altexps_ERCC_percent"))
colSums(as.matrix(filt_reasons))
# low_lib_size            low_n_features high_altexps_ERCC_percent 
# 2                         4                         7 
# discard 
# 9

# filter
storm_sce.filt <- storm_sce.raw[,!filt_reasons$discard]

saveRDS(storm_sce.filt,
        file = "/Users/ben.johnson/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_filt_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")
#storm_sce.filt <- readRDS("/Volumes/projects/shen/projects/2021_02_25_ben_projects/laptop_backup/kb_python/kb_python_storm_test/storm_100k_sce_filt_hg38_ens101.rds")

# how many known enhancers per cell
storm_sce.filt$NumEnhExpressed <- colSums2(counts(altExp(storm_sce.filt, "Enhancers")) > 0)
# re-calculate num genes
storm_sce.filt$NumGenesExpressed <- colSums2(counts(storm_sce.filt) > 0)
# re-calculate depth
storm_sce.filt$depth_size <- colSums2(counts(storm_sce.filt))

library(ggplot2)
library(viridis)

ggplot(colData(storm_sce.filt), aes(x = NumEnhExpressed, y = depth_size, color = depth_size)) +
  geom_point() +
  scale_color_viridis_c(end = 0.8, option = "magma") +
  theme_bw(12)

## subsample to 150k as well to compare
# do it again with all PINTS gro/pro-cap data for K562, distal
storm_counts_files <- list.files("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/150k",
                                 pattern = "UniqueAndMult-EM.mtx",
                                 recursive = TRUE,
                                 full.names = TRUE)
storm_counts_files <- dirname(storm_counts_files)

first_gsub <- gsub("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/storm/150k/",
                   "", storm_counts_files)
names(storm_counts_files) <- gsub("_solo_pintsoutput/GeneFull/raw",
                                  "", first_gsub)

# read in data
storm_counts_raw <- lapply(names(storm_counts_files), function(x) {
  # get the path
  kb_mm <- read_starsolo_gene(storm_counts_files[x],
                              x)
  return(kb_mm)
})

# merge
storm_counts_raw.mat <- do.call(cbind, storm_counts_raw)

# make the SCE object
storm_sce.raw <- SingleCellExperiment(assays = list(counts = storm_counts_raw.mat))

# get the detected genes
# should round these...
storm_sce.raw$NumGenesExpressed <- colSums2(counts(storm_sce.raw) > 0)

# pull out the ERCCs
ercc.sce <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "ERCC") <- ercc.sce
storm_sce.raw <- storm_sce.raw[grep("ERCC-", rownames(storm_sce.raw), invert = TRUE),]

# create a separate altexp for the enhancers
enh.sce <- storm_sce.raw[grep("pints_k562_distal_enh", rownames(storm_sce.raw)),]
altExp(storm_sce.raw, "Enhancers") <- enh.sce
storm_sce.raw <- storm_sce.raw[grep("pints_k562_distal_enh", rownames(storm_sce.raw), invert = TRUE),]

# storm_sce.raw.150k <- readRDS("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_150k_raw_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")

# QC
library(scuttle)
storm_sce.raw.150k <- addPerCellQCMetrics(storm_sce.raw.150k)
filt_reasons <- perCellQCFilters(storm_sce.raw.150k,
                                 sub.fields = c("altexps_ERCC_percent"))
colSums(as.matrix(filt_reasons))
# low_lib_size            low_n_features high_altexps_ERCC_percent                   discard 
# 3                         9                         6                        12 

# filter
storm_sce.filt.150k <- storm_sce.raw.150k[,!filt_reasons$discard]

saveRDS(storm_sce.filt.150k,
        file = "/Users/ben.johnson/Documents/manuscripts/storm_seq/erna/k562_storm/storm/storm_150k_filt_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")

# how many known enhancers per cell
storm_sce.filt.150k$NumEnhExpressed <- colSums2(counts(altExp(storm_sce.filt.150k, "Enhancers")) > 0)
# re-calculate num genes
storm_sce.filt.150k$NumGenesExpressed <- colSums2(counts(storm_sce.filt.150k) > 0)
# re-calculate depth
storm_sce.filt.150k$depth_size <- colSums2(counts(storm_sce.filt.150k))

library(ggplot2)
library(viridis)

ggplot(colData(storm_sce.filt.150k), aes(x = NumEnhExpressed, y = depth_size, color = depth_size)) +
  geom_point() +
  scale_color_viridis_c(end = 0.8, option = "magma") +
  theme_bw(12)

## vasa
## accidentally forgot to change the output folder name from last time...
## so the current output has the pints enhancers now for vasa
## need to change this. current date: 2024-12-11
vasa_counts_files <- list.files("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/vasa_plate/full_depth/150k_to_keep",
                                pattern = "UniqueAndMult-EM.mtx",
                                recursive = TRUE,
                                full.names = TRUE)
vasa_counts_files <- dirname(vasa_counts_files)

first_gsub <- gsub("/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/erna/vasa_plate/full_depth/150k_to_keep/",
                   "", vasa_counts_files)
names(vasa_counts_files) <- gsub("output/GeneFull/raw",
                                 "", first_gsub)

# read in data
vasa_counts_raw <- lapply(names(vasa_counts_files), function(x) {
  # get the path
  kb_mm <- read_starsolo_gene(vasa_counts_files[x],
                              x)
  return(kb_mm)
})

# merge
vasa_counts_raw.mat <- do.call(cbind, vasa_counts_raw)

# make the SCE object
vasa_sce.raw <- SingleCellExperiment(assays = list(counts = vasa_counts_raw.mat))

# get the detected genes
vasa_sce.raw$NumGenesExpressed <- colSums2(counts(vasa_sce.raw) > 0)

# pull out the ERCCs
ercc.sce <- vasa_sce.raw[grep("ERCC-", rownames(vasa_sce.raw)),]
altExp(vasa_sce.raw, "ERCC") <- ercc.sce
vasa_sce.raw <- vasa_sce.raw[grep("ERCC-", rownames(vasa_sce.raw), invert = TRUE),]

# create a separate altexp for the enhancers
enh.sce <- vasa_sce.raw[grep("pints_k562_distal_enh", rownames(vasa_sce.raw)),]
altExp(vasa_sce.raw, "Enhancers") <- enh.sce
vasa_sce.raw <- vasa_sce.raw[grep("pints_k562_distal_enh", rownames(vasa_sce.raw), invert = TRUE),]

# vasa_sce.raw <- readRDS("~/Documents/manuscripts/storm_seq/erna/k562_vasa/vasa_plate_raw_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")

# QC
vasa_sce.raw <- addPerCellQCMetrics(vasa_sce.raw)
filt_reasons <- perCellQCFilters(vasa_sce.raw)
colSums(as.matrix(filt_reasons))
# low_lib_size low_n_features        discard 
# 3              2              3 

# filter
vasa_sce.filt <- vasa_sce.raw[,!filt_reasons$discard]

# remove controls too
p1_neg_controls <- c("vai_vasa_plate1_merged_009", "vai_vasa_plate1_merged_047",
                     "vai_vasa_plate1_merged_108", "vai_vasa_plate1_merged_123",
                     "vai_vasa_plate1_merged_207", "vai_vasa_plate1_merged_318",
                     "vai_vasa_plate1_merged_334", "vai_vasa_plate1_merged_379")
p2_neg_controls <- c("vai_vasa_plate2_merged_009", "vai_vasa_plate2_merged_047",
                     "vai_vasa_plate2_merged_108", "vai_vasa_plate2_merged_123",
                     "vai_vasa_plate2_merged_207", "vai_vasa_plate2_merged_318",
                     "vai_vasa_plate2_merged_334", "vai_vasa_plate2_merged_379")
vasa_sce.filt <- vasa_sce.filt[,!colnames(vasa_sce.filt) %in% c(p1_neg_controls,
                                                                p2_neg_controls)]

saveRDS(vasa_sce.filt,
        file = "/Users/ben.johnson/Documents/manuscripts/storm_seq/erna/k562_vasa/vasa_plate_filt_sce_with_pints_distal_groprocap_k562_enh_annots_hg38_ens101.rds")
#vasa_sce.filt <- readRDS("/Volumes/projects/shen/projects/2021_02_25_ben_projects/laptop_backup/kb_python/kb_python_vasa_test/vasa_100k_sce_filt_hg38_ens101.rds")

# how many known enhancers per cell
vasa_sce.filt$NumEnhExpressed <- colSums2(counts(altExp(vasa_sce.filt, "Enhancers")) > 0)
# re-calculate num genes
vasa_sce.filt$NumGenesExpressed <- colSums2(counts(vasa_sce.filt) > 0)
# re-calculate depth
vasa_sce.filt$depth_size <- colSums2(counts(vasa_sce.filt))

library(ggplot2)
library(viridis)

ggplot(colData(vasa_sce.filt), aes(x = NumEnhExpressed, y = depth_size, color = depth_size)) +
  geom_point() +
  scale_color_viridis_c(end = 0.8, option = "magma") +
  theme_bw(12)

# grab the 150k subsampled bulk k562 gene counts too
k562_rep1_counts <- read.delim("K5621_bulk_star_pints_enhReadsPerGene.out.tab",
                               header = FALSE)
k562_rep1_counts <- as.data.frame(k562_rep1_counts[-c(1:4), c(1,4)])
colnames(k562_rep1_counts) <- c("ens_id", "K562_1")
# count number of detected enhancers
k562_rep1_counts.enh <- k562_rep1_counts[grep("pints", k562_rep1_counts$ens_id),]
table(k562_rep1_counts.enh$K562_1 > 0)
# FALSE  TRUE
# 30124  1542

k562_rep2_counts <- read.delim("K5622_bulk_star_pints_enhReadsPerGene.out.tab",
                               header = FALSE)
k562_rep2_counts <- as.data.frame(k562_rep2_counts[-c(1:4), c(1,4)])
colnames(k562_rep2_counts) <- c("ens_id", "K562_2")
# count number of detected enhancers
k562_rep2_counts.enh <- k562_rep2_counts[grep("pints", k562_rep2_counts$ens_id),]
table(k562_rep2_counts.enh$K562_2 > 0)
# FALSE  TRUE
# 30263  1403

mean(c(1542, 1403))
# 1472.5

bulk_enh <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/bulk/subsample_150k/bulk_k562_reps1_2_pints_distal_enh_counts.txt")
bulk_enh_all <- rowSums(bulk_enh[,c(2:3)])
bulk_enh.filt <- bulk_enh[bulk_enh_all > 0,]

bulk_enh.filt.nonzero <- apply(bulk_enh.filt[,2:3], 1, function(row) all(row != 0))
bulk_enh.filt.ovlp <- bulk_enh.filt[bulk_enh.filt.nonzero,]

# how many unique enhancers per technology?
# storm
table(rowSums2(counts(altExp(storm_sce.filt.150k, "Enhancers"))) > 0)
# FALSE  TRUE 
# 26510  5156

# vasa
table(rowSums2(counts(altExp(vasa_sce.filt, "Enhancers"))) > 0)
# FALSE  TRUE 
# 25868  5798
# whoa... that's an insane complexity...
# how many cells share the same sets of enhancers?
storm_overlapping_enh <- apply(as.matrix(counts(altExp(storm_sce.filt.150k, "Enhancers"))), 1, function(row) sum(row != 0))
storm_overlapping_enh <- storm_overlapping_enh[storm_overlapping_enh > 0]
summary(storm_overlapping_enh)

# randomly subsample the vasa data to the same number of cells as storm
vasa_sce.filt.subsamp <- vasa_sce.filt[,sample(colnames(vasa_sce.filt), size = 100, replace = FALSE)]
vasa_overlapping_enh <- apply(as.matrix(counts(altExp(vasa_sce.filt.subsamp, "Enhancers"))), 1, function(row) sum(row != 0))
vasa_overlapping_enh <- vasa_overlapping_enh[vasa_overlapping_enh > 0]
summary(vasa_overlapping_enh)
# ope. that explains things. most are only found in 1 cell...

storm_overlapping_enh.filt <- storm_overlapping_enh[storm_overlapping_enh > 2]
vasa_overlapping_enh.filt <- vasa_overlapping_enh[vasa_overlapping_enh > 2]

# okay, how many are also found in the bulk?
table(names(storm_overlapping_enh.filt) %in% bulk_enh.filt$ens_id)
# FALSE  TRUE 
# 1554   237 

table(names(vasa_overlapping_enh.filt) %in% bulk_enh.filt$ens_id)
# FALSE  TRUE 
# 631   187

# let's make a complexity curve!
storm_enh_pints <- data.frame(detected = storm_sce.filt$NumEnhExpressed,
                              depth = storm_sce.filt$depth_size,
                              technology = "STORM-seq")
storm_enh_pints.150k <- data.frame(detected = storm_sce.filt.150k$NumEnhExpressed,
                                   depth = storm_sce.filt.150k$depth_size,
                                   technology = "STORM-seq")
vasa_enh_pints <- data.frame(detected = vasa_sce.filt$NumEnhExpressed,
                             depth = vasa_sce.filt$depth_size,
                             technology = "VASA-seq")
all_enh <- rbind(storm_enh_pints.150k,
                 vasa_enh_pints)

# Calculate medians per technology
enh_medians <- all_enh %>%
  group_by(technology) %>%
  summarise(median_detected = median(detected), median_depth = median(depth), .groups = "drop")

ggplot(all_enh, aes(x = detected, y = depth,
                    color = technology, group = technology)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(data = enh_medians, aes(y = median_depth, x = median_detected), 
             color = "white", size = 4, shape = 21, stroke = 1, fill = "black") +
  scale_y_continuous(limits = c(0, 1e5),
                     breaks = seq(0, 1e5, 2.5e4)) +
  scale_x_continuous(limits = c(0, 400),
                     breaks = seq(0, 400, 100)) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF")) +
  ylab("Unique UMI Fragments") +
  xlab("Distal K-562 Enhancers Detected\n(150k reads/cell)") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

# plot out the number of cells where consistent enh
# are detected?
storm_overlapping_enh <- apply(as.matrix(counts(altExp(storm_sce.filt.150k, "Enhancers"))), 1, function(row) sum(row != 0))
vasa_overlapping_enh <- apply(as.matrix(counts(altExp(vasa_sce.filt.subsamp, "Enhancers"))), 1, function(row) sum(row != 0))

enh_ovlp <- data.frame(ovlp = c(storm_overlapping_enh,
                                vasa_overlapping_enh),
                       technology = c(rep("STORM-seq", length(storm_overlapping_enh)),
                                      rep("VASA-seq", length(vasa_overlapping_enh))))
enh_ovlp <- enh_ovlp[enh_ovlp$ovlp > 1,]

# Calculate quantile breakpoints
breaks <- c(2, 5, 10, 15, Inf)

# Bin the data into deciles
binned_data <- data.frame(
  technology = enh_ovlp$technology,
  bin = cut(enh_ovlp$ovlp, breaks = breaks, include.lowest = TRUE)
)

# Summarize counts for each bin and technology
binned_summary <- binned_data %>%
  group_by(technology, bin) %>%
  summarise(count = n(), .groups = "drop")

ggplot(binned_summary, aes(x = bin, y = count,
                           color = technology, group = technology)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF")) +
  scale_x_discrete(labels = c(`[2,5]` = "2-5",
                              `(5,10]` = "5-10",
                              `(10,15]` = "10-15",
                              `(15,Inf]` = ">15")) +
  scale_y_continuous(limits = c(0, 1800),
                     breaks = seq(0, 1800, 600)) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  xlab("Shared Enhancer Expression\n(Cell Count, 150k reads/cell)") +
  ylab("Distal K-562 Enhancers") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

## plot the profiles of STORM, VASA, pro-cap, and TT-seq
# storm
# Load necessary libraries
library(tidyverse)

storm_profile_data.fwd <- read.table("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/merged_storm_k562_solo_pints_aligned.sorted.bpm_norm.fwd.data_matrix.tab",
                                     skip = 3, header = FALSE, sep = "\t")
colnames(storm_profile_data.fwd) <- paste0("Bin_", 1:ncol(storm_profile_data.fwd))
storm_profile_data.rev <- read.table("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/merged_storm_k562_solo_pints_aligned.sorted.bpm_norm.rev.data_matrix.tab",
                                     skip = 3, header = FALSE, sep = "\t")
colnames(storm_profile_data.rev) <- paste0("Bin_", 1:ncol(storm_profile_data.rev))

# filter away all zeroes
storm_profile_data.fwd <- storm_profile_data.fwd[rowSums(storm_profile_data.fwd) > 0,]
storm_profile_data.rev <- storm_profile_data.rev[rowSums(storm_profile_data.rev) > 0,]

# Reshape the data into a long format for ggplot2
storm_profile_data_long.fwd <- storm_profile_data.fwd %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin)) # Ensure Bin is treated as numeric for plotting

storm_profile_data_long.rev <- storm_profile_data.rev %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin))

# Calculate the mean signal for each bin
storm_bin_means.fwd <- storm_profile_data_long.fwd %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))
storm_bin_means.rev <- storm_profile_data_long.rev %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))

# rescale for plotting
storm_bin_means.rev$Mean_Signal.rescale <- storm_bin_means.rev$Mean_Signal/3

# Plot the smoothed average profile using ggplot2
ggplot(storm_bin_means.fwd, aes(x = Bin, y = Mean_Signal)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

ggplot(storm_bin_means.rev, aes(x = Bin, y = -Mean_Signal)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

## vasa
vasa_profile_data.fwd <- read.table("~/Documents/manuscripts/storm_seq/erna/k562_vasa/merged_vasa_starsolo_pints_groprocap_distal_enh.ens101.fwd.pints_groprocap_distal_enh.bpm_norm.data_matrix.tab",
                                     skip = 3, header = FALSE, sep = "\t")
colnames(vasa_profile_data.fwd) <- paste0("Bin_", 1:ncol(vasa_profile_data.fwd))
vasa_profile_data.rev <- read.table("~/Documents/manuscripts/storm_seq/erna/k562_vasa/merged_vasa_starsolo_pints_groprocap_distal_enh.ens101.rev.pints_groprocap_distal_enh.bpm_norm.data_matrix.tab",
                                     skip = 3, header = FALSE, sep = "\t")
colnames(vasa_profile_data.rev) <- paste0("Bin_", 1:ncol(vasa_profile_data.rev))

# filter away all zeroes
vasa_profile_data.fwd <- vasa_profile_data.fwd[rowSums(vasa_profile_data.fwd) > 0,]
vasa_profile_data.rev <- vasa_profile_data.rev[rowSums(vasa_profile_data.rev) > 0,]

# Reshape the data into a long format for ggplot2
vasa_profile_data_long.fwd <- vasa_profile_data.fwd %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin)) # Ensure Bin is treated as numeric for plotting

vasa_profile_data_long.rev <- vasa_profile_data.rev %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin))

# Calculate the mean signal for each bin
vasa_bin_means.fwd <- vasa_profile_data_long.fwd %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))
vasa_bin_means.rev <- vasa_profile_data_long.rev %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))

# rescale for plotting
vasa_bin_means.fwd$Mean_Signal.rescale <- vasa_bin_means.fwd$Mean_Signal/3

# Plot the smoothed average profile using ggplot2
ggplot(vasa_bin_means.fwd, aes(x = Bin, y = Mean_Signal)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

ggplot(vasa_bin_means.rev, aes(x = Bin, y = -Mean_Signal)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

## tt-seq
ttseq_rep1.fwd <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/GSM4610686_L_K562_Rep1.coverage.track.plus_no_chr.pints_distal_grocap_intergenic.data_matrix.tab",
                             skip = 3, sep = '\t', header = FALSE)
ttseq_rep1.rev <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/GSM4610686_L_K562_Rep1.coverage.track.minus_no_chr.pints_distal_grocap_intergenic.data_matrix.tab",
                             skip = 3, sep = '\t', header = FALSE)


ttseq_rep2.fwd <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/GSM4610687_L_K562_Rep2.coverage.track.plus_no_chr.pints_distal_grocap_intergenic.data_matrix.tab",
                             skip = 3, sep = '\t', header = FALSE)
ttseq_rep2.rev <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/GSM4610687_L_K562_Rep2.coverage.track.minus_no_chr.pints_distal_grocap_intergenic.data_matrix.tab",
                             skip = 3, sep = '\t', header = FALSE)

# combine and average
ttseq_all.fwd <- (ttseq_rep1.fwd + ttseq_rep2.fwd)/2
ttseq_all.rev <- (ttseq_rep1.rev + ttseq_rep2.rev)/2

ttseq_all.fwd <- ttseq_all.fwd[rowSums(ttseq_all.fwd) > 0,]
colnames(ttseq_all.fwd) <- paste0("Bin_", 1:ncol(ttseq_all.fwd))

ttseq_all.rev <- ttseq_all.rev[rowSums(ttseq_all.rev) > 0,]
colnames(ttseq_all.rev) <- paste0("Bin_", 1:ncol(ttseq_all.rev))

# Reshape the data into a long format for ggplot2
ttseq_profile_data_long.fwd <- ttseq_all.fwd %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin)) # Ensure Bin is treated as numeric for plotting

ttseq_profile_data_long.rev <- ttseq_all.rev %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin))

# Calculate the mean signal for each bin
ttseq_bin_means.fwd <- ttseq_profile_data_long.fwd %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))
ttseq_bin_means.rev <- ttseq_profile_data_long.rev %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))

# divide signal to put on appropriate scale with BPM norm'd
# STORM data and pro-cap since it's really the shape we are after
ttseq_bin_means.fwd$Mean_Signal.rescale <- ttseq_bin_means.fwd$Mean_Signal/200
ttseq_bin_means.rev$Mean_Signal.rescale <- ttseq_bin_means.rev$Mean_Signal/200

# Plot the smoothed average profile using ggplot2
ggplot(ttseq_bin_means.fwd, aes(x = Bin, y = Mean_Signal.rescale)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )
ggplot(ttseq_bin_means.rev, aes(x = Bin, y = -Mean_Signal.rescale)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

## pro-cap data
procap_plus <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/procap_encode/pro_cap_encode_reps1_2_ave_plus_strand.pints_distal_1kbfilt_expand.data_matrix.tab",
                          skip = 3, sep = '\t', header = FALSE)
procap_minus <- read.delim("~/Documents/manuscripts/storm_seq/erna/k562_storm/procap_encode/pro_cap_encode_reps1_2_ave_minus_strand.pints_distal_1kbfilt_expand.data_matrix.tab",
                           skip = 3, sep = '\t', header = FALSE)

procap_plus <- procap_plus[rowSums(procap_plus) > 0,]
colnames(procap_plus) <- paste0("Bin_", 1:ncol(procap_plus))

procap_minus <- procap_minus[!rowSums(procap_minus) == 0,]
colnames(procap_minus) <- paste0("Bin_", 1:ncol(procap_minus))

# Reshape the data into a long format for ggplot2
procap_plus_profile_data_long <- procap_plus %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin)) # Ensure Bin is treated as numeric for plotting

procap_minus_profile_data_long <- procap_minus %>%
  pivot_longer(cols = starts_with("Bin_"), 
               names_to = "Bin", 
               names_prefix = "Bin_", 
               values_to = "Signal") %>%
  mutate(Bin = as.numeric(Bin)) # Ensure Bin is treated as numeric for plotting

# Calculate the mean signal for each bin
procap_plus_bin_means <- procap_plus_profile_data_long %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))

procap_minus_bin_means <- procap_minus_profile_data_long %>%
  group_by(Bin) %>%
  summarise(Mean_Signal = mean(Signal, na.rm = TRUE))

# rescale for plotting
procap_plus_bin_means$Mean_Signal.rescale <- procap_plus_bin_means$Mean_Signal * 6
procap_minus_bin_means$Mean_Signal.rescale <- procap_minus_bin_means$Mean_Signal * 6

# Plot the smoothed average profile using ggplot2
ggplot(procap_plus_bin_means, aes(x = Bin, y = Mean_Signal.rescale)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

ggplot(procap_minus_bin_means, aes(x = Bin, y = Mean_Signal.rescale)) +
  geom_line(color = "blue", size = 1) +  # Line through the mean
  # geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +  # Optional smoothing
  theme_minimal() +
  labs(
    title = "Smoothed Average Signal Profile",
    x = "Bins",
    y = "Mean Signal",
    caption = "Generated from deeptools computeMatrix output"
  )

## intersect the storm and vasa data?
# NOTE: the forward vasa data is the reverse storm and vice versa!
# this is due to the strand differences between the techs
vasa_profile_data.fwd.storm <- vasa_profile_data.fwd[rownames(vasa_profile_data.fwd) %in% rownames(storm_profile_data.rev),]
vasa_profile_data.rev.storm <- vasa_profile_data.rev[rownames(vasa_profile_data.rev) %in% rownames(storm_profile_data.fwd),]

## combine and plot
to_plot <- data.frame(bin = c(storm_bin_means.fwd$Bin,
                              storm_bin_means.rev$Bin,
                              vasa_bin_means.rev$Bin,
                              vasa_bin_means.fwd$Bin,
                              ttseq_bin_means.fwd$Bin,
                              ttseq_bin_means.rev$Bin,
                              procap_plus_bin_means$Bin,
                              procap_minus_bin_means$Bin),
                      mean_signal = c(storm_bin_means.fwd$Mean_Signal - min(storm_bin_means.fwd$Mean_Signal) + 0.005,
                                      -(storm_bin_means.rev$Mean_Signal.rescale - min(storm_bin_means.rev$Mean_Signal.rescale) + 0.005),
                                      vasa_bin_means.rev$Mean_Signal - min(vasa_bin_means.rev$Mean_Signal) + 0.005,
                                      -(vasa_bin_means.fwd$Mean_Signal.rescale - min(vasa_bin_means.fwd$Mean_Signal.rescale) + 0.005),
                                      ttseq_bin_means.fwd$Mean_Signal.rescale - min(ttseq_bin_means.fwd$Mean_Signal.rescale) + 0.005,
                                      -(ttseq_bin_means.rev$Mean_Signal.rescale - min(ttseq_bin_means.rev$Mean_Signal.rescale) + 0.005),
                                      procap_plus_bin_means$Mean_Signal.rescale - min(procap_plus_bin_means$Mean_Signal.rescale) + 0.005,
                                      procap_minus_bin_means$Mean_Signal.rescale - max(procap_minus_bin_means$Mean_Signal.rescale) - 0.005),
                      assay = c(rep("STORM-seq +", length(storm_bin_means.fwd$Bin)),
                                rep("STORM-seq -", length(storm_bin_means.rev$Bin)),
                                rep("VASA-seq +", length(vasa_bin_means.rev$Bin)),
                                rep("VASA-seq -", length(vasa_bin_means.fwd$Bin)),
                                rep("TT-seq +", length(ttseq_bin_means.fwd$Bin)),
                                rep("TT-seq -", length(ttseq_bin_means.rev$Bin)),
                                rep("Pro-cap +", length(procap_plus_bin_means$Bin)),
                                rep("Pro-cap -", length(procap_minus_bin_means$Bin))))

to_plot$assay <- factor(to_plot$assay, levels = c("STORM-seq +",
                                                  "STORM-seq -",
                                                  "VASA-seq +",
                                                  "VASA-seq -",
                                                  "TT-seq +",
                                                  "TT-seq -",
                                                  "Pro-cap +",
                                                  "Pro-cap -"))
ggplot(to_plot, aes(x = bin, y = mean_signal,
                    group = assay, color = assay,
                    alpha = assay, linewidth = assay)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(`STORM-seq +`="#63197FFF",
                                `STORM-seq -`="#63197FFF",
                                `VASA-seq +`="#CF5917FF",
                                `VASA-seq -`="#CF5917FF",
                                `TT-seq +` = "#2E6E8EFF",
                                `TT-seq -` = "#2E6E8EFF",
                                `Pro-cap +` = "#000004FF",
                                `Pro-cap -` = "#000004FF")) +
  scale_linewidth_manual(values = c(`STORM-seq +` = 1.4,
                                    `STORM-seq -` = 1.4,
                                    `VASA-seq +` = 1.4,
                                    `VASA-seq -` = 1.4,
                                    `TT-seq +` = 0.9,
                                    `TT-seq -` = 0.9,
                                    `Pro-cap +` = 0.9,
                                    `Pro-cap -` = 0.9)) +
  scale_alpha_manual(values = c(`STORM-seq +` = 1,
                                `STORM-seq -` = 1,
                                `VASA-seq +` = 1,
                                `VASA-seq -` = 1,
                                `TT-seq +` = 0.5,
                                `TT-seq -` = 0.5,
                                `Pro-cap +` = 0.5,
                                `Pro-cap -` = 0.5)) +
  scale_x_continuous(labels = c(`0` = "-0.5kb",
                              #`25` = "-0.25kb",
                              `50` = "Center",
                              #`75` = "+0.25kb",
                              `100` = "+0.5kb"),
                     breaks = seq(0,100,50)) +
  # scale_y_continuous(limits = c(-0.05, 0.05),
  #                    breaks = seq(-0.05, 0.05, 0.025)) +
  geom_vline(xintercept = 50, linewidth = 1,
             linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 1,
             linetype = "solid") +
  ylab("Normalized Coverage") +
  theme_bw(12) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(size = 18)
  )

# number unique enhancers covered by each tech
length(union(rownames(storm_profile_data.fwd), rownames(storm_profile_data.rev)))
# 3767
length(union(rownames(ttseq_all.fwd), rownames(ttseq_all.rev)))
# 8088
length(union(rownames(procap_plus), rownames(procap_minus)))
# 7599
length(union(rownames(vasa_profile_data.fwd), rownames(vasa_profile_data.rev)))

# might also be useful to show capture of stable/unstable
# transcript distributions?
