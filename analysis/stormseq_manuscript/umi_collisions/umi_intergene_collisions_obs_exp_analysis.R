## pull in the real collision rates for obs/exp calcs
# NOTE: it doesn't matter what is used for the 16 bp UMI as it is
## effectively 0 across the board
## NOTE: use the 8k as the exp for storm, vasa, and ss3xpress
storm_collisions.hek <- read.delim("inter_gene_collisions_observed/storm/storm_100k_umi_collision_rate_hek.txt",
                               header = FALSE)
ss3_collisions <- read.delim("inter_gene_collisions_observed/ss3x/ss3xpress_100k_umi_collision_rates.txt",
                             header = FALSE)
vasa_collisions.hek <- read.delim("inter_gene_collisions_observed/vasa/vasa_100k_umi_collision_rate_hek.txt",
                              header = FALSE)
sstotal_collisions <- read.delim("inter_gene_collisions_observed/smart_seq_total/sstotal_umi_collisions.txt",
                                 header = FALSE)

# upfront filter and recalc
# remove low quality sstotal cells based on fragments
sstotal_collisions <- sstotal_collisions[sstotal_collisions$V3 > 4000,]

# recalculate vasa rate since we only have 4096 with 6 x N UMI/UFI
vasa_collisions.hek$recalc <- 1.0 - vasa_collisions.hek$V4

# re-ingest the simulation rates for the exp rate at 8k
umi_6_res.long <- readRDS("inter_gene_collision_simulations/simulation_results/umi_6/umi_6_simulation_results.rds")
umi_8_res.long <- readRDS("inter_gene_collision_simulations/simulation_results/umi_8/umi_8_simulation_results.rds")
umi_10_res.long <- readRDS("inter_gene_collision_simulations/simulation_results/umi_10/umi_10_simulation_results.rds")
umi_16_res.long <- readRDS("inter_gene_collision_simulations/simulation_results/umi_16/umi_16_simulation_results.rds")

# melt
umi_6_res.long <- reshape2::melt(umi_6_res.long, id.vars = rownames(umi_6_res.long))
umi_8_res.long <- reshape2::melt(umi_8_res.long, id.vars = rownames(umi_8_res.long))
umi_10_res.long <- reshape2::melt(umi_10_res.long, id.vars = rownames(umi_10_res.long))
umi_16_res.long <- reshape2::melt(umi_16_res.long, id.vars = rownames(umi_16_res.long))

# annotate 
umi_6_res.long$umi <- "NNNNNN"
umi_8_res.long$umi <- "NNNNNNNN"
umi_10_res.long$umi <- "NNNNNNNNWW"
umi_16_res.long$umi <- "NNNNNNNNNNNNNNNN"

## obs/exp
mean_8k_6bp <- umi_6_res.long[umi_6_res.long$Var1 %in% "8000",]
mean_8k_6bp <- mean(mean_8k_6bp$value)

mean_8k_8bp <- umi_8_res.long[umi_8_res.long$Var1 %in% "8000",]
mean_8k_8bp <- mean(mean_8k_8bp$value)

mean_8k_10bp <- umi_10_res.long[umi_10_res.long$Var1 %in% "8000",]
mean_8k_10bp <- mean(mean_8k_10bp$value)

mean_8k_16bp <- umi_16_res.long[umi_16_res.long$Var1 %in% "8000",]
# mean_8k_16bp <- mean(mean_8k_16bp) # the exp rate is 0.0 at 16 x N UMI
mean_8k_16bp <- 0.0

storm_collisions.hek$obs_exp <- storm_collisions.hek$V4/mean_8k_8bp
storm_collisions.hek$tech <- "STORM-seq"
vasa_collisions.hek$obs_exp <- vasa_collisions.hek$recalc/mean_8k_6bp
vasa_collisions.hek$tech <- "VASA-seq"
# add a very small pseudocount
sstotal_collisions$obs_exp <- (sstotal_collisions$V4 + 0.001)/(mean_8k_16bp + 0.001)
sstotal_collisions$tech <- "Smart-seq-total"
ss3_collisions$obs_exp <- ss3_collisions$V4/mean_8k_10bp
ss3_collisions$tech <- "Smart-seq3xpress"

# collate
all_collisions <- rbind(storm_collisions.hek[,5:6],
                        vasa_collisions.hek[,6:7],
                        sstotal_collisions[,5:6],
                        ss3_collisions[,5:6])

all_collisions$tech <- factor(all_collisions$tech,
                              levels = c("STORM-seq",
                                         "VASA-seq",
                                         "Smart-seq3xpress",
                                         "Smart-seq-total"))

# plot
library(ggplot2)
# truncate the axis and will remove a few large smart-seq-total points
ggplot(all_collisions, aes(x = tech, y = obs_exp, color = tech)) +
  geom_jitter(width = 0.1, height = 0.1) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF",
                                `Smart-seq-total` = "#8FC2FD",
                                `Smart-seq3xpress` = "#8CCC98"),
                     name = "Technology") +
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymax = ..y.., ymin = ..y.., group = tech),
    width = 0.2,
    position = position_dodge(width = 0.75),
    size = 2,
    color = "black"
  ) +
  ylim(c(0, 10)) +
  geom_hline(yintercept = 1,
             linetype = "dotted",
             color = "black") +
  ylab("Obs/Exp Inter-gene UMI Collision Rate\n(100k reads/cell)") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )
