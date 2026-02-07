## gene fusion analysis
## all done at 150k reads/cell or sample in the case of bulk
## all settings were identical for STAR-Fusion invocation
## followed the tutorial: https://github.com/STAR-Fusion/STAR-Fusion/wiki/STAR-Fusion-scRNA-seq

# known gene fusions from CCLE
k562_fusions <- read.csv("~/Documents/manuscripts/storm_seq/star_fusion/K562 fusions.csv")
k562_fusions <- unique(k562_fusions$Fusion.Name)

# use the calls from the very same sample and the above as ground truth
k562_mrna_fusions <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/bulk/mrna/k562_mrna_star-fusion.fusion_predictions.abridged.tsv")
k562_mrna_fusions.filt <- k562_mrna_fusions[k562_mrna_fusions$X.FusionName %in% k562_fusions,]

# # grab cosmic fusions
# cosmic_samples <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/Cosmic_Sample_Tsv_v100_GRCh38/Cosmic_Sample_v100_GRCh38.tsv.gz")
# cosmic_samples.k562 <- cosmic_samples[cosmic_samples$SAMPLE_NAME %in% c("K-562", "K562"),]
# # somehow there are 54 (!!) unique COSMIC sample IDs for K562...
# cosmic_fusions <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/Cosmic_Fusion_Tsv_v100_GRCh38/Cosmic_Fusion_v100_GRCh38.tsv.gz")
# cosmic_fusions.k562 <- cosmic_fusions[cosmic_fusions$COSMIC_SAMPLE_ID %in% cosmic_samples.k562$COSMIC_SAMPLE_ID,]
# # well this is odd... there aren't _any_ gene fusions in cosmic for K562

# bulk
# full depth
bulk_starF_fd <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/bulk/total/bulk_k562_full_depth_star_fusion_agg_deconv.fusions.abridged.tsv")
bulk_starF_fd.filt <- bulk_starF_fd[bulk_starF_fd$X.FusionName %in% k562_fusions,]
#bulk_1 <- bulk_starF_fd[1:17,]
# there's a header at index 18
#bulk_2 <- bulk_starF_fd[19:33,]
#bulk_starF_fd.intersect <- intersect(bulk_1$X.FusionName, bulk_2$X.FusionName)

# validate the presence of these other fusions in depmap
# thanks chatGPT for the next couple lines
# Split each string into two parts using '--' as the delimiter
# split_vector <- strsplit(bulk_starF_fd.intersect, "--")
# Reverse the order of the parts and concatenate them back together
# reversed_vector <- sapply(split_vector, function(x) paste(rev(x), collapse="--"))

# depmap fusions
# depmap_fusions <- read.csv("~/Documents/manuscripts/storm_seq/star_fusion/depmap_24q2_rel_OmicsFusionFiltered.csv")

# bulk_starF_fd.intersect.depmap <- bulk_starF_fd.intersect[bulk_starF_fd.intersect %in% depmap_fusions$FusionName]
# bulk_starF_fd.intersect.rev.depmap <- reversed_vector[reversed_vector %in% depmap_fusions$FusionName]

# manually add the fusions that are found in depmap
# bulk_starF_fd.curated <- c(bulk_starF_fd.intersect.depmap,
#                            "KANSL1--ARL17A")
# k562_fusions.curated <- union(k562_fusions, bulk_starF_fd.curated)

#150k
# bulk_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/bulk/bulk_k562_150k_star_fusion_agg_deconv.fusions.abridged.tsv")
# bulk_starF.filt <- bulk_starF[bulk_starF$X.FusionName %in% k562_fusions.curated,]

# storm
storm_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/storm/storm_star_fusion_150k_depth_k562_agg_deconv.fusions.abridged.tsv")
storm_starF.filt <- storm_starF[storm_starF$X.FusionName %in% k562_fusions,]

# vasa
vasa_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/vasa_plate/star_fusion_vasa_plateSE_150k_agg_deconv.fusions.abridged.tsv")
vasa_starF.filt <- vasa_starF[vasa_starF$X.FusionName %in% k562_fusions,]

# ss3xpress
ss3_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/ss3/ss3xpress_internal_150k_depth_star_fusion_agg_deconv.fusions.abridged.tsv")
ss3_starF.filt <- ss3_starF[ss3_starF$X.FusionName %in% k562_fusions,]

# calculate number of cells with presence of known fusions
# bulk
num_cells.bulk.fd <- 2
# num_cells.bulk <- length(unique(bulk_starF.filt$Cell))

# storm
num_cells.storm <- length(unique(storm_starF.filt$Cell))
num_cells.storm <- num_cells.storm/112

# vasa
num_cells.vasa <- length(unique(vasa_starF.filt$Cell))
num_cells.vasa <- num_cells.vasa/188

# ss3xpress
num_cells.ss3 <- length(unique(ss3_starF.filt$Cell))
num_cells.ss3 <- num_cells.ss3/692

# calculate the number of fusions per cell we can get
# getNumUniqFusions <- function(df, known.fusions = NULL) {
#   # extract unique fusions
#   sv.uniq <- unique(df[,1])
#   if (!is.null(known.fusions)) {
#     num_known <- length(sv.uniq[sv.uniq %in% known.fusions])
#     return(data.frame(num_known = num_known,
#                       num_unknown = length(sv.uniq)-num_known))
#   } else {
#     return(data.frame(num_uniq = length(sv.uniq)))
#   }
# }
# 
# storm_starF.uniq_fusions <- getNumUniqFusions(storm_starF,
#                                               known.fusions = k562_fusions)
# bulk_starF_fd.uniq_fusions <- getNumUniqFusions(bulk_starF_fd,
#                                                 known.fusions = k562_fusions)
# ss3_starF.uniq_fusions <- getNumUniqFusions(ss3_starF,
#                                             known.fusions = k562_fusions)

# storm
storm_percell_svs <- unlist(lapply(k562_fusions, function(x) {
  res <- storm_starF[storm_starF$X.FusionName %in% x,]
  res.dedup <- res[!duplicated(res$Cell),]
  return(nrow(res.dedup))
}))
names(storm_percell_svs) <- k562_fusions
storm_percell_svs <- (storm_percell_svs/112) * 100

# vasa
vasa_percell_svs <- unlist(lapply(k562_fusions, function(x) {
  res <- vasa_starF[vasa_starF$X.FusionName %in% x,]
  res.dedup <- res[!duplicated(res$Cell),]
  return(nrow(res.dedup))
}))
names(vasa_percell_svs) <- k562_fusions
vasa_percell_svs <- (vasa_percell_svs/188) * 100

# ss3xpress
ss3_percell_svs <- unlist(lapply(k562_fusions, function(x) {
  res <- ss3_starF[ss3_starF$X.FusionName %in% x,]
  res.dedup <- res[!duplicated(res$Cell),]
  return(nrow(res.dedup))
}))
names(ss3_percell_svs) <- k562_fusions
ss3_percell_svs <- (ss3_percell_svs/692) * 100

# gather and weight proportions by FFPM within bulk total and mRNA (CCLE)
library(dplyr)
k562_fusions.ref <- read.csv("~/Documents/manuscripts/storm_seq/star_fusion/K562 fusions.csv")
k562_ref_fusion_ffpm <- k562_fusions.ref %>%
  group_by(Fusion.Name) %>%
  summarise(mean_FFPM = mean(Ffpm, na.rm = TRUE))
k562_ref_fusion_ffpm$method <- "mRNA"

k562_mrna_fusion_ffpm <- k562_mrna_fusions.filt %>%
  group_by(X.FusionName) %>%
  summarise(mean_FFPM = mean(FFPM, na.rm = TRUE))
k562_mrna_fusion_ffpm$method <- "mRNA"

bulk_starF_fd.filt$FFPM <- as.numeric(bulk_starF_fd.filt$FFPM)
k562_total_fusion_ffpm <- bulk_starF_fd.filt %>%
  group_by(X.FusionName) %>%
  summarise(mean_FFPM = mean(FFPM, na.rm = TRUE))
k562_total_fusion_ffpm$method <- "total"

all(k562_mrna_fusion_ffpm$X.FusionName == k562_total_fusion_ffpm$X.FusionName)
# TRUE
ffpm_ratios <- k562_mrna_fusion_ffpm$mean_FFPM / k562_total_fusion_ffpm$mean_FFPM
names(ffpm_ratios) <- k562_total_fusion_ffpm$X.FusionName

# BAG6--SLC44A4      BCR--ABL1 C16orf87--ORC6  IMMP2L--DOCK4   NUP214--XKR3   UPF3A--CDC16    XACT--LRCH2 
# 5.9934126      1.0154589      2.8472222      1.3144126      2.6089983      3.0028925      0.9030149

# weight the per cell svs by FFPM within relative bulk method
# total
# storm
storm_percell_svs.filt <- storm_percell_svs[names(storm_percell_svs) %in% k562_total_fusion_ffpm$X.FusionName]
storm_percell_svs_weighted <- unlist(lapply(names(storm_percell_svs.filt), function(x) {
  prop <- storm_percell_svs.filt[x]
  weight <- ffpm_ratios[x]
  return(prop * weight)
}))

vasa_percell_svs.filt <- vasa_percell_svs[names(vasa_percell_svs) %in% k562_total_fusion_ffpm$X.FusionName]
vasa_percell_svs_weighted <- unlist(lapply(names(vasa_percell_svs.filt), function(x) {
  prop <- vasa_percell_svs.filt[x]
  weight <- ffpm_ratios[x]
  return(prop * weight)
}))

ss3_percell_svs.filt <- ss3_percell_svs[names(ss3_percell_svs) %in% k562_mrna_fusion_ffpm$X.FusionName]
ss3_percell_svs_weighted <- unlist(lapply(names(ss3_percell_svs.filt), function(x) {
  prop <- ss3_percell_svs.filt[x]
  weight <- ffpm_ratios[x]
  return(prop / weight)
}))

# do a bar plot of proportion of cells with known fusions
library(ggplot2)
library(viridis)

all_sv_fracs <- data.frame(technology = c(rep("STORM-seq", 7),
                                          rep("VASA-seq", 7),
                                          rep("Smart-seq3xpress", 7)),
                           fusion_fracs = c(storm_percell_svs_weighted,
                                            vasa_percell_svs_weighted,
                                            ss3_percell_svs_weighted),
                           fusions = names(storm_percell_svs_weighted))

# filter out things that are zeroes
# all_sv_fracs <- all_sv_fracs[all_sv_fracs$fusion_fracs > 0,]
all_sv_fracs$fusions <- gsub("--", "-", all_sv_fracs$fusions)

# remove AC006453.1-FAM230C since it's found in 1 cell in ss3xpress
# all_sv_fracs <- all_sv_fracs[!all_sv_fracs$fusions %in% "AC006453.1-FAM230C",]

# order by storm data
storm_data <- all_sv_fracs[all_sv_fracs$technology == "STORM-seq", ]
ordered_fusions <- storm_data$fusions[order(-storm_data$fusion_fracs)]
all_sv_fracs$fusions <- factor(all_sv_fracs$fusions, levels = ordered_fusions)
all_sv_fracs$technology <- factor(all_sv_fracs$technology,
                                  levels = c("STORM-seq",
                                             "VASA-seq",
                                             "Smart-seq3xpress"))

ggplot(all_sv_fracs, aes (x = technology, y = fusion_fracs, fill = technology)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,20)) +
  scale_fill_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF",
                                `Smart-seq3xpress` = "#8CCC98"),
                     name = "Technology") +
  facet_wrap(~ fusions) +
  ylab("Percent Cells with Fusion\n(FFPM weighted, 150k reads/cell)") +
  theme_bw(12) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "None",
    # legend.title = element_blank(),
    # legend.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14)
  )

# circos plot time...
library(circlize)

## turn the above into a function...
extractSVCoords <- function(cell_to_extract = NULL, sv_table,
                            padding = 10, known_fusions = NULL) {
  ## extract the cell or not
  if (!is.null(cell_to_extract)) {
    sv_table.sub <- sv_table[sv_table$Cell == cell_to_extract,]
  } else {
    sv_table.sub <- sv_table
  }
  
  ## generate the beds
  cell_bed1 <- data.frame(chr = sapply(strsplit(sv_table.sub$LeftBreakpoint, ":"), `[`, 1),
                          start = as.integer(sapply(strsplit(sv_table.sub$LeftBreakpoint, ":"), `[`, 2)),
                          end = as(as.integer(sapply(strsplit(sv_table.sub$LeftBreakpoint, ":"), `[`, 2)) + padding, "integer"),
                          value1 = as(sv_table.sub$JunctionReadCount, "integer"),
                          value2 = as(sv_table.sub$SpanningFragCount, "integer"),
                          fusion = sv_table.sub$X.FusionName)
  cell_bed2 <- data.frame(chr = sapply(strsplit(sv_table.sub$RightBreakpoint, ":"), `[`, 1),
                          start = as.integer(sapply(strsplit(sv_table.sub$RightBreakpoint, ":"), `[`, 2)),
                          end = as(as.integer(sapply(strsplit(sv_table.sub$RightBreakpoint, ":"), `[`, 2)) + padding, "integer"),
                          value1 = as(sv_table.sub$JunctionReadCount, "integer"),
                          value2 = as(sv_table.sub$SpanningFragCount, "integer"),
                          fusion = sv_table.sub$X.FusionName)
  cell_bed1 <- subset(cell_bed1, !duplicated(cell_bed1$fusion))
  cell_bed2 <- subset(cell_bed2, !duplicated(cell_bed2$fusion))
  
  if (!is.null(known_fusions)) {
    # subset to just the known vector of fusions
    cell_bed1 <- cell_bed1[cell_bed1$fusion %in% known_fusions,]
    cell_bed2 <- cell_bed2[cell_bed2$fusion %in% known_fusions,]
  }
  
  return(list(bed1 = cell_bed1,
              bed2 = cell_bed2))
}

# subset to the known 7 fusions we have listed above
# storm
storm_sv_coords <- extractSVCoords(sv_table = storm_starF.filt,
                                   known_fusions = k562_total_fusion_ffpm$X.FusionName)

# initialize some colors
library(viridis)
cols <- viridis_pal(option = "magma", end = 0.8)(7)
names(cols) <- c("BCR--ABL1", "NUP214--XKR3",
                 "BAG6--SLC44A4", "XACT--LRCH2",
                 "C16orf87--ORC6", "UPF3A--CDC16",
                 "IMMP2L--DOCK4")
cols.m <- match(storm_sv_coords[[1]]$fusion,
                names(cols))
my_cols <- cols[cols.m]

circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr6", "chr7", "chr9",
                                                   "chr13", "chr16", "chr22",
                                                   "chrX"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(storm_sv_coords[[1]], storm_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% storm_sv_coords[[1]]$fusion])
circos.clear()

# ss3
ss3_sv_coords <- extractSVCoords(sv_table = ss3_starF.filt,
                                   known_fusions = k562_total_fusion_ffpm$X.FusionName)
# NOTE: need to re-order the fusions to match the colors above.
# they are identical in the number of fusions detected so it's not really a problem
# to re-use the storm circos and the bulk total circos
# make it explicitly clear though
ss3_sv_coords.m <- match(storm_sv_coords[[1]]$fusion,
                         ss3_sv_coords[[1]]$fusion)
ss3_sv_coords <- lapply(ss3_sv_coords, function(x) {
  return(x[ss3_sv_coords.m,])
})
circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr6", "chr7", "chr9",
                                                   "chr13", "chr16", "chr22",
                                                   "chrX"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(ss3_sv_coords[[1]], ss3_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% ss3_sv_coords[[1]]$fusion])
circos.clear()

# vasa
vasa_sv_coords <- extractSVCoords(sv_table = vasa_starF.filt,
                                   known_fusions = k562_total_fusion_ffpm$X.FusionName)
vasa_sv_coords.m <- match(storm_sv_coords[[1]]$fusion,
                         vasa_sv_coords[[1]]$fusion)
vasa_sv_coords <- lapply(vasa_sv_coords, function(x) {
  return(x[vasa_sv_coords.m,])
})
circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr6", "chr7", "chr9",
                                                   "chr13", "chr16", "chr22",
                                                   "chrX"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(vasa_sv_coords[[1]], vasa_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% vasa_sv_coords[[1]]$fusion])
circos.clear()

# bulk total
bulk_total_sv_coords <- extractSVCoords(sv_table = bulk_starF_fd.filt,
                                        known_fusions = k562_total_fusion_ffpm$X.FusionName)
bulk_total_sv_coords.m <- match(storm_sv_coords[[1]]$fusion,
                          bulk_total_sv_coords[[1]]$fusion)
bulk_total_sv_coords <- lapply(bulk_total_sv_coords, function(x) {
  return(x[bulk_total_sv_coords.m,])
})
circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr6", "chr7", "chr9",
                                                   "chr13", "chr16", "chr22",
                                                   "chrX"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(bulk_total_sv_coords[[1]], bulk_total_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% bulk_total_sv_coords[[1]]$fusion])
circos.clear()

## process the full depth RMG-2 bulk mRNA RNA and STORM
## note - though "total" appears in the variable name, it is mrna
rmg2_bulk_total_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/bulk/mrna/rmg2_mrna_star-fusion.fusion_predictions.abridged.tsv")

# filter away the gtex recurrent calls
rmg2_bulk_total_starF.filt <- rmg2_bulk_total_starF[grep("GTEx_recurrent_StarF2019",
                                                         rmg2_bulk_total_starF$annots,
                                                         invert = TRUE),]

# filter away things that are inter- or intrachromosomal that does _not_ have large
# anchor support or is unknown to existing gene fusion databases
rmg2_high_conf_fusions <- rmg2_bulk_total_starF.filt[rmg2_bulk_total_starF.filt$LargeAnchorSupport %in% "YES_LDAS",]
rmg2_potential_fusions <- rmg2_bulk_total_starF.filt[rmg2_bulk_total_starF.filt$LargeAnchorSupport %in% "NO_LDAS",]
# filter to known fusions
starF_fusion_annot.lib <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/fusion_annot_lib.gz",
                                     header = FALSE)
rmg2_potential_fusions.known <- rmg2_potential_fusions[rmg2_potential_fusions$X.FusionName %in% starF_fusion_annot.lib$V1,]
# filter to fusions in known oncogenes
rmg2_potential_fusions.onco <- rmg2_potential_fusions[grep("Oncogene",
                                                           rmg2_potential_fusions$annots),]
# combine 
rmg2_curated_fusions <- rbind(rmg2_high_conf_fusions,
                              rmg2_potential_fusions.known,
                              rmg2_potential_fusions.onco)

# extract the detected fusions
rmg2_fusions <- unique(rmg2_curated_fusions$X.FusionName)

# storm full depth
# storm
rmg2_storm_starF <- read.delim("~/Documents/manuscripts/storm_seq/star_fusion/storm/storm_rmg2_star_fusion_agg_deconv_full_depth.fusions.abridged.tsv")
rmg2_storm_starF.filt <- rmg2_storm_starF[rmg2_storm_starF$X.FusionName %in% unique(rmg2_bulk_total_starF.filt$X.FusionName),]

# pull in the filtered rmg2 object to see how many cells we have
storm_sce.filt.100k <- readRDS("~/Documents/manuscripts/storm_seq/sensitivity_calc/storm_100k_sce_filt_hg38_ens101.rds")
storm_rmg2 <- colnames(storm_sce.filt.100k)[storm_sce.filt.100k$cell_type %in% "RMG2"]
# note - the fusion calls are the _old_ well names and need to be remapped
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
                          header = FALSE)
storm_remap.rmg2.m <- match(storm_rmg2,
                            storm_remap$V2)
storm_remap.rmg2 <- storm_remap[storm_remap.rmg2.m,]
table(unique(rmg2_storm_starF.filt$Cell) %in% storm_remap.rmg2$V1)

# filter to passing cells with the previous QC (e.g. those used in figure 2)
rmg2_storm_starF.filt.cellfilt <- rmg2_storm_starF.filt[rmg2_storm_starF.filt$Cell %in% storm_remap.rmg2$V1,]

rmg2_storm_percell_svs <- unlist(lapply(rmg2_fusions, function(x) {
  res <- rmg2_storm_starF.filt.cellfilt[rmg2_storm_starF.filt.cellfilt$X.FusionName %in% x,]
  res.dedup <- res[!duplicated(res$Cell),]
  return(nrow(res.dedup))
}))
names(rmg2_storm_percell_svs) <- rmg2_fusions
# LAMA5--GATA5      CLTC--FP671120.4            ITCH--ASIP           TET3--DGUOK        CDH1--CYP4F24P 
# 68                     1                    61                    18                     3 
# AC093821.1--LINC01091        STK33--DENND2B           LIPC--MTMR3 
# 7                             12                        15
rmg2_storm_percell_svs <- (rmg2_storm_percell_svs/97) * 100
# LAMA5--GATA5      CLTC--FP671120.4            ITCH--ASIP           TET3--DGUOK        CDH1--CYP4F24P 
# 70.103093              1.030928             62.886598             18.556701              3.092784 
# AC093821.1--LINC01091        STK33--DENND2B           LIPC--MTMR3 
# 7.216495                      12.371134             15.463918

# plot out per cell fractions by fusion for STORM-seq
names(rmg2_storm_percell_svs) <- gsub("--", "-", names(rmg2_storm_percell_svs))

# make it a df
rmg2_sv_fracs <- data.frame(technology = rep("STORM-seq", length(rmg2_storm_percell_svs)),
                            fusion_fracs = rmg2_storm_percell_svs,
                            fusions = names(rmg2_storm_percell_svs))
rmg2_sv_fracs <- rmg2_sv_fracs[order(rmg2_sv_fracs$fusion_fracs, decreasing = TRUE),]
rmg2_sv_fracs$fusions <- factor(rmg2_sv_fracs$fusions,
                                levels = rmg2_sv_fracs$fusions)
ggplot(rmg2_sv_fracs, aes (x = fusions, y = fusion_fracs, fill = technology)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,20)) +
  scale_fill_manual(values = c(`STORM-seq`="#63197FFF"),
                    name = "Technology") +
  ylab("Percent Cells with Fusion") +
  ggtitle("Per-cell RMG-II Gene Fusions in STORM-seq") +
  theme_bw(12) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 90,
                               vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "None"
    # legend.title = element_blank(),
    # legend.text = element_text(size = 16),
    # strip.background = element_rect(fill = "white", color = "black"),
    # strip.text = element_text(size = 14)
  )

# circos time
storm_rmg2_sv_coords <- extractSVCoords(sv_table = rmg2_storm_starF.filt.cellfilt,
                                        known_fusions = rmg2_fusions)
# initialize some colors
library(viridis)
cols <- viridis_pal(option = "viridis", end = 0.8)(8)
names(cols) <- names(rmg2_storm_percell_svs)
cols.m <- match(storm_rmg2_sv_coords[[1]]$fusion,
                names(cols))
my_cols <- cols[cols.m]

circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr2", "chr4", "chr11",
                                                   "chr15", "chr16", "chr17",
                                                   "chr19", "chr20", "chr21",
                                                   "chr22"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(storm_rmg2_sv_coords[[1]], storm_rmg2_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% storm_rmg2_sv_coords[[1]]$fusion])
circos.clear()

# bulk rna
bulk_rmg2_sv_coords <- extractSVCoords(sv_table = rmg2_bulk_total_starF.filt,
                                        known_fusions = rmg2_fusions)
bulk_total_sv_coords.m <- match(storm_rmg2_sv_coords[[1]]$fusion,
                                bulk_rmg2_sv_coords[[1]]$fusion)
bulk_rmg2_sv_coords <- lapply(bulk_rmg2_sv_coords, function(x) {
  return(x[bulk_total_sv_coords.m,])
})

circos.initializeWithIdeogram(species = "hg38",
                              chromosome.index = c("chr2", "chr4", "chr11",
                                                   "chr15", "chr16", "chr17",
                                                   "chr19", "chr20", "chr21",
                                                   "chr22"),
                              labels.cex = 2, axis.labels.cex = 1,
                              major.by = 5e7, sort.chr = TRUE)
circos.genomicLink(bulk_rmg2_sv_coords[[1]], bulk_rmg2_sv_coords[[2]],
                   lwd = 3, col = my_cols[names(my_cols) %in% bulk_rmg2_sv_coords[[1]]$fusion])
circos.clear()
