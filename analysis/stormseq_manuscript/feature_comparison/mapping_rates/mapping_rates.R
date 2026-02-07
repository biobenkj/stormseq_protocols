## compare mapping rates in HEK cells across technologies

#STORM
storm_mapping_rate <- read.delim("~/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/mapping_rate_storm_with_filenames.txt",
                                 header = FALSE)
storm_mapping_rate <- storm_mapping_rate[1:384,]
storm_mapping_rate$V1 <- gsub("_kb_tcc", "", storm_mapping_rate$V1)
storm_mapping_rate$platform <- "STORM-seq"

# need to remap the data to only get HEK cells
# remap the wells to the proper place and remove controls
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
                          header = FALSE)
storm_remap.match <- match(storm_mapping_rate$V1,
                           storm_remap$V1)
storm_remap <- storm_remap[storm_remap.match,]

#careful
all(storm_remap$V1 == storm_mapping_rate$V1)
# TRUE

# re-map the names/wells
storm_mapping_rate$V1 <- storm_remap$V2

# only keep HEK cells
hek293t_cells <- c(paste0("B", 1:24),
                   paste0("E", 1:24),
                   paste0("G", 1:6),
                   paste0("G", 8:24),
                   paste0("I", 1:6),
                   paste0("I", 8:24),
                   paste0("L", c(1,3:24)),
                   paste0("N", 1:24))
names(hek293t_cells) <- rep("HEK293T", length(hek293t_cells))
cell_type.match <- match(storm_mapping_rate$V1,
                         hek293t_cells)
storm_mapping_rate.hek <- storm_mapping_rate[storm_mapping_rate$V1 %in% hek293t_cells,]

# only keep cells that pass QC filters
storm_sce.filt <- readRDS("~/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/storm_100k_sce_filt_hg38_ens101.rds")
storm_mapping_rate.hek <- storm_mapping_rate.hek[storm_mapping_rate.hek$V1 %in% colnames(storm_sce.filt),]


#SMART-seq Total
sstotal_mapping_rate <- read.delim("~/Documents/manuscripts/storm_seq/subsample_analysis/sstotal_hek/mapping_rate_sstotal_with_filenames.txt",
                                 header = FALSE)
sstotal_mapping_rate$platform <- "SMART-seq Total"
sstotal_mapping_rate$V1 <- gsub("_kb_tcc", "", sstotal_mapping_rate$V1)
# only keep cells that pass QC filters
sstotal_sce.filt <- readRDS("~/Documents/manuscripts/storm_seq/subsample_analysis/sstotal_hek/sstotal_100k_sce_filt_hg38_ens101.rds")

sstotal_mapping_rate <- sstotal_mapping_rate[sstotal_mapping_rate$V1 %in% colnames(sstotal_sce.filt),]

#VASA
vasa_mapping_rate <- read.delim("~/Documents/manuscripts/storm_seq/subsample_analysis/vasa_hek/mapping_rate_vasa_with_filenames.txt",
                                   header = FALSE)
p1_neg_controls <- c("VAI-MR-v001_HNVCGBGXN_S3_009", "VAI-MR-v001_HNVCGBGXN_S3_047",
                     "VAI-MR-v001_HNVCGBGXN_S3_108", "VAI-MR-v001_HNVCGBGXN_S3_123",
                     "VAI-MR-v001_HNVCGBGXN_S3_207", "VAI-MR-v001_HNVCGBGXN_S3_318",
                     "VAI-MR-v001_HNVCGBGXN_S3_334", "VAI-MR-v001_HNVCGBGXN_S3_379")
p2_neg_controls <- c("VAI-MR-v002_HNVCGBGXN_S4_009", "VAI-MR-v002_HNVCGBGXN_S4_047",
                     "VAI-MR-v002_HNVCGBGXN_S4_108", "VAI-MR-v002_HNVCGBGXN_S4_123",
                     "VAI-MR-v002_HNVCGBGXN_S4_207", "VAI-MR-v002_HNVCGBGXN_S4_318",
                     "VAI-MR-v002_HNVCGBGXN_S4_334", "VAI-MR-v002_HNVCGBGXN_S4_379")

vasa_cell_type <- c(rep("K562", 192),
                    rep("HEK293T", 192),
                    rep("K562", 192),
                    rep("HEK293T", 192))
vasa_mapping_rate$cell_type <- vasa_cell_type
vasa_mapping_rate.hek <- vasa_mapping_rate[vasa_mapping_rate$cell_type %in% "HEK293T",]
# remove any controls
vasa_mapping_rate.hek <- vasa_mapping_rate.hek[!gsub("_kb_tcc", "", vasa_mapping_rate.hek$V1) %in% c(p1_neg_controls,
                                                                                                     p2_neg_controls),]

vasa_mapping_rate.hek <- vasa_mapping_rate.hek[,c(1:2)]
vasa_mapping_rate.hek$platform <- "VASA-seq"

# only keep cells passing QC
vasa_sce.filt <- readRDS("~/Documents/manuscripts/storm_seq/subsample_analysis/vasa_hek/vasa_100k_sce_filt_hg38_ens101.rds")

vasa_mapping_rate.hek <- vasa_mapping_rate.hek[gsub("_kb_tcc", "", vasa_mapping_rate.hek$V1) %in% colnames(vasa_sce.filt),]

library(ggplot2)
library(ggridges)

all_mapping_rates <- rbind(storm_mapping_rate.hek,
                           sstotal_mapping_rate,
                           vasa_mapping_rate.hek)
#all_mapping_rates$V2 <- as.numeric(all_mapping_rates$V2)
all_mapping_rates$platform <- factor(all_mapping_rates$platform,
                                     levels = c("STORM-seq",
                                                "SMART-seq Total",
                                                "VASA-seq"))

ggplot(all_mapping_rates, aes(x = V2, y = platform, fill = platform)) +
  geom_density_ridges() +
  scale_fill_manual(values = c(`STORM-seq` = "#641A80FF",
                               `SMART-seq Total` = "#8FC2FDFF",
                               `VASA-seq` = "#CF5A17")) +
  xlim(c(0,100)) +
  theme_bw(12) +
  ggtitle("Mapping rates across single cells") +
  theme(
    legend.position = "None",
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              size = 18),
    axis.text = element_text(size = 18,
                             color = "black")
  )
