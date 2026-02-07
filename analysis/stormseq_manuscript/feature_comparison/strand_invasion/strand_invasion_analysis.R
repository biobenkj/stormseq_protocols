# calculate strand invasion per technology using the techniques
# described in smart-seq3xpress paper
# code was directly provided by Michael and Rickard
# converted to python with help from chatGPT

## load libraries
library(dplyr)
library(ggplot2)

## storm
storm_si <- read.delim("~/Documents/manuscripts/storm_seq/strand_invasion/storm/sandberg_r2py_umi_results_storm_updated.txt")

# filter away anything that wasn't put in the UB tag by star solo
storm_si.umi_filt <- storm_si[!storm_si$UB %in% "-",]

# remap things to pull out the hek cells only for comparison
storm_remap <- read.delim("/Users/ben.johnson/Documents/manuscripts/storm_seq/subsample_analysis/storm_hek/well_map.txt",
                          header = FALSE)

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

k562_cells <- c(paste0("A", 1:24),
                paste0("D", 1:24),
                paste0("F", 1:24),
                paste0("K", 1:24),
                paste0("P", 1:24))
names(k562_cells) <- rep("K562", length(k562_cells))

rmg2_cells <- c(paste0("C", 1:24),
                paste0("H", c(1,3:24)),
                paste0("J", 1:24),
                paste0("M", 1:24),
                paste0("O", 1:24))
names(rmg2_cells) <- rep("RMG2", length(rmg2_cells))

# add a column for just the well name to remap
storm_si.umi_filt$well_id <- gsub("_soloAligned.sortedByCoord.out",
                                  "",
                                  storm_si.umi_filt$RG)

storm_remap.m <- match(storm_si.umi_filt$well_id,
                       storm_remap$V1)
storm_remap <- storm_remap[storm_remap.m,]
all(storm_remap$V1 == storm_si.umi_filt$well_id)
#TRUE
storm_si.umi_filt$remap <- storm_remap$V2
storm_si.umi_filt.hek <- storm_si.umi_filt[storm_si.umi_filt$remap %in% hek293t_cells,]
storm_si.umi_filt.k562 <- storm_si.umi_filt[storm_si.umi_filt$remap %in% k562_cells,]
storm_si.umi_filt.rmg2 <- storm_si.umi_filt[storm_si.umi_filt$remap %in% rmg2_cells,]

# summarize
storm_si.umi_filt.hek.summary <- storm_si.umi_filt.hek %>%
  group_by(remap) %>%
  summarise(
    total_entries = n(),
    match_check_true = sum(match_check == "True"),
    match_check_1mm_true = sum(match_check_1mm == "True"),
    match_check_2mm_true = sum(match_check_2mm == "True"),
    match_check_gg_true = sum(match_check_gg == "True")
  ) %>%
  mutate(
    match_check_prop = (match_check_true / total_entries) * 100,
    match_check_1mm_prop = (match_check_1mm_true / total_entries) * 100,
    match_check_2mm_prop = (match_check_2mm_true / total_entries) * 100,
    match_check_gg_prop = (match_check_gg_true / total_entries) * 100
  )

storm_si.umi_filt.k562.summary <- storm_si.umi_filt.k562 %>%
  group_by(remap) %>%
  summarise(
    total_entries = n(),
    match_check_true = sum(match_check == "True"),
    match_check_1mm_true = sum(match_check_1mm == "True"),
    match_check_2mm_true = sum(match_check_2mm == "True"),
    match_check_gg_true = sum(match_check_gg == "True")
  ) %>%
  mutate(
    match_check_prop = (match_check_true / total_entries) * 100,
    match_check_1mm_prop = (match_check_1mm_true / total_entries) * 100,
    match_check_2mm_prop = (match_check_2mm_true / total_entries) * 100,
    match_check_gg_prop = (match_check_gg_true / total_entries) * 100
  )

storm_si.umi_filt.rmg2.summary <- storm_si.umi_filt.rmg2 %>%
  group_by(remap) %>%
  summarise(
    total_entries = n(),
    match_check_true = sum(match_check == "True"),
    match_check_1mm_true = sum(match_check_1mm == "True"),
    match_check_2mm_true = sum(match_check_2mm == "True"),
    match_check_gg_true = sum(match_check_gg == "True")
  ) %>%
  mutate(
    match_check_prop = (match_check_true / total_entries) * 100,
    match_check_1mm_prop = (match_check_1mm_true / total_entries) * 100,
    match_check_2mm_prop = (match_check_2mm_true / total_entries) * 100,
    match_check_gg_prop = (match_check_gg_true / total_entries) * 100
  )

# toss out the low quality cells
storm_si.umi_filt.hek.summary <- storm_si.umi_filt.hek.summary[storm_si.umi_filt.hek.summary$total_entries > 1000,]
summary(storm_si.umi_filt.hek.summary$match_check_1mm_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.860   2.394   2.558   2.565   2.717   3.969 

sd(storm_si.umi_filt.hek.summary$match_check_1mm_prop)
# 0.2768306
mean(storm_si.umi_filt.hek.summary$match_check_1mm_prop)
# 2.565122

storm_si.umi_filt.k562.summary <- storm_si.umi_filt.k562.summary[storm_si.umi_filt.k562.summary$total_entries > 10000,]
summary(storm_si.umi_filt.k562.summary$match_check_1mm_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.686   2.651   2.894   2.862   3.094   3.462 

storm_si.umi_filt.rmg2.summary <- storm_si.umi_filt.rmg2.summary[storm_si.umi_filt.rmg2.summary$total_entries > 10000,]
summary(storm_si.umi_filt.rmg2.summary$match_check_1mm_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.815   2.552   2.751   2.757   2.988   3.463 

# ss3xpress
ss3_si <- read.delim("~/Documents/manuscripts/storm_seq/strand_invasion/ss3/sandberg_r2py_umi_results_hek_100k_updated.txt")

# filter away anything that wasn't put in the UB tag by star solo
ss3_si.umi_filt <- ss3_si[!ss3_si$UB %in% "-",]

# add in a well ID/cell barcode
ss3_si.umi_filt$cbc <- gsub("Aligned.sortedByCoord.out",
                            "",
                            ss3_si.umi_filt$RG)

# summarize
ss3_si.umi_filt.hek.summary <- ss3_si.umi_filt %>%
  group_by(cbc) %>%
  summarise(
    total_entries = n(),
    match_check_true = sum(match_check == "True"),
    match_check_1mm_true = sum(match_check_1mm == "True"),
    match_check_2mm_true = sum(match_check_2mm == "True"),
    match_check_gg_true = sum(match_check_gg == "True")
  ) %>%
  mutate(
    match_check_prop = (match_check_true / total_entries) * 100,
    match_check_1mm_prop = (match_check_1mm_true / total_entries) * 100,
    match_check_2mm_prop = (match_check_2mm_true / total_entries) * 100,
    match_check_gg_prop = (match_check_gg_true / total_entries) * 100
  )

summary(ss3_si.umi_filt.hek.summary$match_check_1mm_prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.387   2.102   2.192   2.133   2.268   2.475

sd(ss3_si.umi_filt.hek.summary$match_check_1mm_prop)
# 0.2468112
mean(ss3_si.umi_filt.hek.summary$match_check_1mm_prop)
# 2.13253

# bring in the K-562 as well but currently can't access HPC...

# NOTE: cannot actually look at strand invasion for VASA and SStotal UMI
# because the reads with the UMI/UFI don't have cDNA to map and check...

# plot out
# compare ss3x and storm
# show across cell types
library(ggplot2)
library(viridis)

ss3_to_plot <- data.frame(zero_mismatch = ss3_si.umi_filt.hek.summary$match_check_prop,
                          one_mismatch = ss3_si.umi_filt.hek.summary$match_check_1mm_prop,
                          technology = rep("Smart-seq3xpress", length(ss3_si.umi_filt.hek.summary$match_check_prop)),
                          cell_type = rep("HEK293T", length(ss3_si.umi_filt.hek.summary$match_check_prop)))
storm_to_plot.hek <- data.frame(zero_mismatch = storm_si.umi_filt.hek.summary$match_check_prop,
                                one_mismatch = storm_si.umi_filt.hek.summary$match_check_1mm_prop,
                                technology = rep("STORM-seq", length(storm_si.umi_filt.hek.summary$match_check_prop)),
                                cell_type = rep("HEK293T", length(storm_si.umi_filt.hek.summary$match_check_prop)))
storm_to_plot.k562 <- data.frame(zero_mismatch = storm_si.umi_filt.k562.summary$match_check_prop,
                                one_mismatch = storm_si.umi_filt.k562.summary$match_check_1mm_prop,
                                technology = rep("STORM-seq", length(storm_si.umi_filt.k562.summary$match_check_prop)),
                                cell_type = rep("K-562", length(storm_si.umi_filt.k562.summary$match_check_prop)))
storm_to_plot.rmg2 <- data.frame(zero_mismatch = storm_si.umi_filt.rmg2.summary$match_check_prop,
                                 one_mismatch = storm_si.umi_filt.rmg2.summary$match_check_1mm_prop,
                                 technology = rep("STORM-seq", length(storm_si.umi_filt.rmg2.summary$match_check_prop)),
                                 cell_type = rep("RMG-2", length(storm_si.umi_filt.rmg2.summary$match_check_prop)))

to_plot <- rbind(ss3_to_plot,
                 storm_to_plot.hek,
                 storm_to_plot.k562,
                 storm_to_plot.rmg2)

to_plot.melt <- reshape2::melt(to_plot, id.vars = c("technology", "cell_type"))

# Calculate medians
median_data <- to_plot.melt %>%
  group_by(technology, cell_type, variable) %>%
  summarize(median_value = median(value), .groups = "drop")

# Plot
ggplot(median_data, aes(x = variable, y = median_value, color = technology, group = technology)) +
  geom_point(size = 3) + # Median points
  geom_line(size = 1) +  # Line connecting median points
  facet_wrap(~ cell_type) + # Facet by cell type
  scale_y_continuous(limits = c(0,16),
                     breaks = seq(0,16,2)) +
  scale_color_manual(values = c("Smart-seq3xpress" = "#8CCC98",
                                "STORM-seq" = "#63197FFF")) +
  labs(
    title = "Median Strand Invasion",
    x = "Variable",
    y = "Percent Strand Invasion",
    color = "Technology"
  ) +
  geom_hline(yintercept = c(15, 2.5), linetype = "dashed",
             color = "black", linewidth = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )


ggplot(to_plot.melt, aes(x = cell_type, y = value,
                         color = technology)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(data = to_plot.medians, aes(y = value, x = variable), 
             color = "white", size = 4, shape = 21, stroke = 1, fill = "black") +
  # scale_y_continuous(limits = c(0, 1e5),
  #                    breaks = seq(0, 1e5, 2.5e4)) +
  # scale_x_continuous(limits = c(0, 400),
  #                    breaks = seq(0, 400, 100)) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `Smart-seq3xpress` = "#8CCC98")) +
  facet_wrap(technology~variable) +
  # ylab("Unique UMI Fragments") +
  # xlab("Distal K-562 Enhancers Detected\n(150k reads/cell)") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )



