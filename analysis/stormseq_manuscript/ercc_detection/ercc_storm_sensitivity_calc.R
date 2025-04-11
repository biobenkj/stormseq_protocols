## calculate detection/capture efficiency across cell types using ERCCs
## ercc info is hard to find
## used data from here: https://jdblischak.github.io/singleCellSeq/analysis/ercc-counts.html
## and here: https://github.com/jdblischak/singleCellSeq/tree/master/data
storm_sce.filt.100k <- readRDS("storm_100k_sce_filt_hg38_ens101.rds")
storm_sce.filt.250k <- readRDS("storm_250k_filt_sce_ens101.rds")
storm_sce.filt.500k <- readRDS("storm_500k_filt_sce_gene_kb_ens101.rds")
storm_sce.filt.1m <- readRDS("storm_1M_filt_sce_ens101_kb.rds")

ercc_info.annot <- read.delim("ercc-info.txt",
                              skip = 1, row.names = 1, header = FALSE)
colnames(ercc_info.annot) <- c("id", "subgroup", "conc_mix1", "conc_mix2",
                               "expected_fc", "log2_mix1_mix2")
# we only have mix 1 here!
# calculate number of molecules per cell at 1:1M dilution
ercc_info.annot <- ercc_info.annot[order(ercc_info.annot$id),]
# according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4191988/
# there are 28000 ERCC molecules per cell if used at 1:20000 dilution
# figure 1d legend and methods for k562 for these numbers
# we have 1:1M so... 
molecules_per_cell_stock <- 28000/50
ratio <- molecules_per_cell_stock / sum(ercc_info.annot[,"conc_mix1"])
ercc_info.annot$counts <- ercc_info.annot[,"conc_mix1"] * ratio
# find ERCCs that are on about 1 molecule per reaction
# this works out to be 1.26 molecules at this dilution.
ercc_info.annot.onemolecule <- ercc_info.annot[ercc_info.annot$counts >=1 & ercc_info.annot$counts <=2,]

# what if we look just below 1?
ercc_info.annot[ercc_info.annot$counts >=0.5 & ercc_info.annot$counts <=1,]
# id subgroup conc_mix1 conc_mix2 expected_fc log2_mix1_mix2    counts
# 29 ERCC-00035        B  117.1875 117.18750        1.00           0.00 0.6339659
# 52 ERCC-00044        C  117.1875 175.78125        0.67          -0.58 0.6339659
# 7  ERCC-00095        A  117.1875  29.29688        4.00           2.00 0.6339659
# 75 ERCC-00112        D  117.1875 234.37500        0.50          -1.00 0.6339659
# 8  ERCC-00131        A  117.1875  29.29688        4.00           2.00 0.6339659

# okay, at most, we should be able to see 21 species of ERCCs at this dilution
# with a minimum of 1.26 molecules/cell 
ercc_info.annot.onemolecule_min <- ercc_info.annot[ercc_info.annot$counts >=1,]

# what is our detection rate of non-zero erccs at 100k?
ercc_storm_100k <- altExp(storm_sce.filt.100k, "ERCC")
ercc_storm_100k$NumGenesExpressed <- colSums2(counts(ercc_storm_100k) > 0)
plot(hist(ercc_storm_100k$NumGenesExpressed))
summary(ercc_storm_100k$NumGenesExpressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 21.00   29.00   31.00   30.63   33.00   39.00 

# well it looks like we are slightly above. how many would we expect if we threshold
# at 0.5 molecules and above?
ercc_info.annot.halfmolecule_min <- ercc_info.annot[ercc_info.annot$counts >=0.5,]
# 26 molecules... that's close to the mean!
# keep in mind we are subsampling here from a high depth to 100k
# let's just see how many we find per cell
ercc_ovlp_100k <- lapply(colnames(ercc_storm_100k), function(x) {
  ercc_sub <- counts(ercc_storm_100k[,x])
  ercc_names_detected <- names(ercc_sub[ercc_sub[,1] > 0,])
  tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
  fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  # calculate correlation with expected counts
  ercc_sub <- ercc_sub[ercc_names_detected,]
  ercc_info.sub <- ercc_info.annot[ercc_info.annot$id %in% ercc_names_detected,]
  ercc_mod <- summary(lm(ercc_sub ~ ercc_info.sub$counts))
  return(data.frame(cell_id = x,
                    ppv = precision,
                    sensitivity = recall,
                    adj.rsq = ercc_mod$r.squared,
                    ercc_sub = ercc_sub,
                    ercc_info_counts = ercc_info.sub$counts))
})
ercc_ovlp_100k <- do.call(rbind, ercc_ovlp_100k)

ercc_ovlp_250k <- lapply(colnames(ercc_storm_250k), function(x) {
  ercc_sub <- counts(ercc_storm_250k[,x])
  ercc_names_detected <- names(ercc_sub[ercc_sub[,1] > 0,])
  tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
  fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  # calculate correlation with expected counts
  ercc_sub <- ercc_sub[ercc_names_detected,]
  ercc_info.sub <- ercc_info.annot[ercc_info.annot$id %in% ercc_names_detected,]
  ercc_mod <- summary(lm(ercc_sub ~ ercc_info.sub$counts))
  return(data.frame(cell_id = x,
                    ppv = precision,
                    sensitivity = recall,
                    adj.rsq = ercc_mod$r.squared,
                    ercc_sub = ercc_sub,
                    ercc_info_counts = ercc_info.sub$counts))
})
ercc_ovlp_250k <- do.call(rbind, ercc_ovlp_250k)

ercc_ovlp_500k <- lapply(colnames(ercc_storm_500k), function(x) {
  ercc_sub <- counts(ercc_storm_500k[,x])
  ercc_names_detected <- names(ercc_sub[ercc_sub[,1] > 0,])
  tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
  fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  # calculate correlation with expected counts
  ercc_sub <- ercc_sub[ercc_names_detected,]
  ercc_info.sub <- ercc_info.annot[ercc_info.annot$id %in% ercc_names_detected,]
  ercc_mod <- summary(lm(ercc_sub ~ ercc_info.sub$counts))
  return(data.frame(cell_id = x,
                    ppv = precision,
                    sensitivity = recall,
                    adj.rsq = ercc_mod$r.squared,
                    ercc_sub = ercc_sub,
                    ercc_info_counts = ercc_info.sub$counts))
})
ercc_ovlp_500k <- do.call(rbind, ercc_ovlp_500k)

ercc_ovlp_1m <- lapply(colnames(ercc_storm_1m), function(x) {
  ercc_sub <- counts(ercc_storm_1m[,x])
  ercc_names_detected <- names(ercc_sub[ercc_sub[,1] > 0,])
  tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
  fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  # calculate correlation with expected counts
  ercc_sub <- ercc_sub[ercc_names_detected,]
  ercc_info.sub <- ercc_info.annot[ercc_info.annot$id %in% ercc_names_detected,]
  ercc_mod <- summary(lm(ercc_sub ~ ercc_info.sub$counts))
  return(data.frame(cell_id = x,
                    ppv = precision,
                    sensitivity = recall,
                    adj.rsq = ercc_mod$r.squared,
                    ercc_sub = ercc_sub,
                    ercc_info_counts = ercc_info.sub$counts))
})
ercc_ovlp_1m <- do.call(rbind, ercc_ovlp_1m)

# example obs/exp for paper
hek_obs_exp_cells <- c("B10", "E14", "B4", "N20")
hek_obs_exp.cells.100k <- ercc_ovlp_100k[ercc_ovlp_100k$cell_id %in% hek_obs_exp_cells,]
hek_obs_exp.cells.250k <- ercc_ovlp_250k[ercc_ovlp_250k$cell_id %in% hek_obs_exp_cells,]
hek_obs_exp.cells.500k <- ercc_ovlp_500k[ercc_ovlp_500k$cell_id %in% hek_obs_exp_cells,]
hek_obs_exp.cells.1m <- ercc_ovlp_1m[ercc_ovlp_1m$cell_id %in% hek_obs_exp_cells,]

hek_obs_exp.cells.100k$ercc_info_counts[hek_obs_exp.cells.100k$ercc_info_counts < 1] <- 1
hek_obs_exp.cells.250k$ercc_info_counts[hek_obs_exp.cells.250k$ercc_info_counts < 1] <- 1
hek_obs_exp.cells.500k$ercc_info_counts[hek_obs_exp.cells.500k$ercc_info_counts < 1] <- 1
hek_obs_exp.cells.1m$ercc_info_counts[hek_obs_exp.cells.1m$ercc_info_counts < 1] <- 1

# annotate
hek_obs_exp.cells.100k$depth <- "100k"
hek_obs_exp.cells.250k$depth <- "250k"
hek_obs_exp.cells.500k$depth <- "500k"
hek_obs_exp.cells.1m$depth <- "1M"

hek_obs_exp.cells <- rbind(hek_obs_exp.cells.100k,
                           hek_obs_exp.cells.250k,
                           hek_obs_exp.cells.500k,
                           hek_obs_exp.cells.1m)

hek_obs_exp.cells.100k.sub <- hek_obs_exp.cells.100k[!hek_obs_exp.cells.100k$cell_id %in% "E14",]

library(ggplot2)
# rotate through each cell as needed as plotting needed to be modified so
# couldn't just facet for the figure.
ggplot(hek_obs_exp.cells.100k.sub, aes(x = log2(ercc_info_counts + 1), y = log2(ercc_sub + 1), color = depth)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  # ylim(c(0,12)) +
  # xlim(c(0,12)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(0,12), breaks = c(seq(0,12,2))) +
  scale_y_continuous(limits = c(0,12), breaks = c(seq(0,12,2))) +
  facet_wrap(~ cell_id, nrow = 1, ncol = 3) +
  ylab("Observed ERCC Counts\n(log2(counts) + 1)") +
  xlab("Expected ERCC Counts (log2(counts) + 1)") +
  scale_color_manual(values = c(`100k` = "black")) +
                                # `250k` = "gold",
                                # `500k` = "red",
                                # `1M` = "black")) +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none"
  )

  

library(dplyr)
# Define manual thresholds
manual_thresholds <- c(0.001, 0.01, 0.1, 1, 5, 10, 25, 50)

# Initialize a list to store results
ercc_manual_results.100k <- list()
ercc_manual_results.250k <- list()
ercc_manual_results.500k <- list()
ercc_manual_results.1m <- list()

# Loop over each threshold to calculate precision and recall
for (threshold in manual_thresholds) {
  ercc_ovlp <- lapply(colnames(ercc_storm_100k), function(x) {
    ercc_sub <- counts(ercc_storm_100k[, x])
    ercc_names_detected <- names(ercc_sub[ercc_sub[, 1] >= threshold, ])  # Apply manual threshold
    
    tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
    fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    
    data.frame(
      cell_id = x,
      detection_threshold = threshold,
      ppv = precision,
      sensitivity = recall
    )
  })
  
  ercc_manual_results.100k[[as.character(threshold)]] <- do.call(rbind, ercc_ovlp)
}

# Combine results into a single data frame
ercc_manual_results_df.100k <- do.call(rbind, ercc_manual_results.100k)

# 250k
for (threshold in manual_thresholds) {
  ercc_ovlp <- lapply(colnames(ercc_storm_250k), function(x) {
    ercc_sub <- counts(ercc_storm_250k[, x])
    ercc_names_detected <- names(ercc_sub[ercc_sub[, 1] >= threshold, ])  # Apply manual threshold
    
    tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
    fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    
    data.frame(
      cell_id = x,
      detection_threshold = threshold,
      ppv = precision,
      sensitivity = recall
    )
  })
  
  ercc_manual_results.250k[[as.character(threshold)]] <- do.call(rbind, ercc_ovlp)
}

# Combine results into a single data frame
ercc_manual_results_df.250k <- do.call(rbind, ercc_manual_results.250k)

# 500k
for (threshold in manual_thresholds) {
  ercc_ovlp <- lapply(colnames(ercc_storm_500k), function(x) {
    ercc_sub <- counts(ercc_storm_500k[, x])
    ercc_names_detected <- names(ercc_sub[ercc_sub[, 1] >= threshold, ])  # Apply manual threshold
    
    tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
    fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    
    data.frame(
      cell_id = x,
      detection_threshold = threshold,
      ppv = precision,
      sensitivity = recall
    )
  })
  
  ercc_manual_results.500k[[as.character(threshold)]] <- do.call(rbind, ercc_ovlp)
}

# Combine results into a single data frame
ercc_manual_results_df.500k <- do.call(rbind, ercc_manual_results.500k)

# 1m
for (threshold in manual_thresholds) {
  ercc_ovlp <- lapply(colnames(ercc_storm_1m), function(x) {
    ercc_sub <- counts(ercc_storm_1m[, x])
    ercc_names_detected <- names(ercc_sub[ercc_sub[, 1] >= threshold, ])  # Apply manual threshold
    
    tp <- sum(ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    fn <- sum(!ercc_info.annot.halfmolecule_min$id %in% ercc_names_detected)
    fp <- sum(!ercc_names_detected %in% ercc_info.annot.halfmolecule_min$id)
    
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
    
    data.frame(
      cell_id = x,
      detection_threshold = threshold,
      ppv = precision,
      sensitivity = recall
    )
  })
  
  ercc_manual_results.1m[[as.character(threshold)]] <- do.call(rbind, ercc_ovlp)
}

# Combine results into a single data frame
ercc_manual_results_df.1m <- do.call(rbind, ercc_manual_results.1m)

# cell type
k562_wells <- rownames(colData(storm_sce.filt.100k)[storm_sce.filt.100k$cell_type %in% "K562",])
hek293t_wells <- rownames(colData(storm_sce.filt.100k)[storm_sce.filt.100k$cell_type %in% "HEK293T",])
rmg2_wells <- rownames(colData(storm_sce.filt.100k)[storm_sce.filt.100k$cell_type %in% "RMG2",])

k562_ercc_ovlp_100k <- ercc_ovlp_100k[ercc_ovlp_100k$cell_id %in% k562_wells,]
k562_ercc_ovlp_100k$cell_type <- "K562"
hek_ercc_ovlp_100k <- ercc_ovlp_100k[ercc_ovlp_100k$cell_id %in% hek293t_wells,]
hek_ercc_ovlp_100k$cell_type <- "HEK293T"
rmg2_ercc_ovlp_100k <- ercc_ovlp_100k[ercc_ovlp_100k$cell_id %in% rmg2_wells,]
rmg2_ercc_ovlp_100k$cell_type <- "RMG2"

k562_ercc_ovlp_250k <- ercc_ovlp_250k[ercc_ovlp_250k$cell_id %in% k562_wells,]
k562_ercc_ovlp_250k$cell_type <- "K562"
hek_ercc_ovlp_250k <- ercc_ovlp_250k[ercc_ovlp_250k$cell_id %in% hek293t_wells,]
hek_ercc_ovlp_250k$cell_type <- "HEK293T"
rmg2_ercc_ovlp_250k <- ercc_ovlp_250k[ercc_ovlp_250k$cell_id %in% rmg2_wells,]
rmg2_ercc_ovlp_250k$cell_type <- "RMG2"

k562_ercc_ovlp_500k <- ercc_ovlp_500k[ercc_ovlp_500k$cell_id %in% k562_wells,]
k562_ercc_ovlp_500k$cell_type <- "K562"
hek_ercc_ovlp_500k <- ercc_ovlp_500k[ercc_ovlp_500k$cell_id %in% hek293t_wells,]
hek_ercc_ovlp_500k$cell_type <- "HEK293T"
rmg2_ercc_ovlp_500k <- ercc_ovlp_500k[ercc_ovlp_500k$cell_id %in% rmg2_wells,]
rmg2_ercc_ovlp_500k$cell_type <- "RMG2"

k562_ercc_ovlp_1m <- ercc_ovlp_1m[ercc_ovlp_1m$cell_id %in% k562_wells,]
k562_ercc_ovlp_1m$cell_type <- "K562"
hek_ercc_ovlp_1m <- ercc_ovlp_1m[ercc_ovlp_1m$cell_id %in% hek293t_wells,]
hek_ercc_ovlp_1m$cell_type <- "HEK293T"
rmg2_ercc_ovlp_1m <- ercc_ovlp_1m[ercc_ovlp_1m$cell_id %in% rmg2_wells,]
rmg2_ercc_ovlp_1m$cell_type <- "RMG2"

# add depth
all_to_plot.100k <- rbind(k562_ercc_ovlp_100k,
                     hek_ercc_ovlp_100k,
                     rmg2_ercc_ovlp_100k)
all_to_plot.100k$depth <- "100k"

all_to_plot.250k <- rbind(k562_ercc_ovlp_250k,
                          hek_ercc_ovlp_250k,
                          rmg2_ercc_ovlp_250k)
all_to_plot.250k$depth <- "250k"

all_to_plot.500k <- rbind(k562_ercc_ovlp_500k,
                          hek_ercc_ovlp_500k,
                          rmg2_ercc_ovlp_500k)
all_to_plot.500k$depth <- "500k"

all_to_plot.1m <- rbind(k562_ercc_ovlp_1m,
                     hek_ercc_ovlp_1m,
                     rmg2_ercc_ovlp_1m)
all_to_plot.1m$depth <- "1M"

# plot things out
library(ggplot2)

all_to_plot <- rbind(all_to_plot.100k,
                     all_to_plot.250k,
                     all_to_plot.500k,
                     all_to_plot.1m)

all_to_plot$depth <- factor(all_to_plot$depth,
                            levels = c("100k",
                                       "250k",
                                       "500k",
                                       "1M"))

# sensitivity first
ggplot(all_to_plot, aes(x = cell_type, y = sensitivity, color = depth)) +
  # geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_viridis_d(option = "inferno", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
  ggtitle("Sensitivity Across Cell Types for STORM-seq",
          subtitle = "Detection of Expected ERCCs at ≥1 Molecule") +
  ylab("Sensitivity") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

# adj r^2 to expected counts
ggplot(all_to_plot, aes(x = cell_type, y = adj.rsq, color = depth)) +
  #geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_viridis_d(option = "inferno", end = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
  ggtitle("Observed vs. Expected ERCC Counts Per Cell",
          subtitle = "Expected ERCCs at ≥1 Molecule") +
  ylab("Adjusted R^2") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )