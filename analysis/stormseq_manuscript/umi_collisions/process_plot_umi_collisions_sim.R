## plot out the inter-gene collision rates

## simulation
umi_lens <- c("umi_10", "umi_16",
              "umi_6", "umi_8")

## gene numbers
genes_sampled <- c("100", "1000",
                   "10000", "12000",
                   "2000", "4000",
                   "500", "6000",
                   "8000")

## loop through and get all 10 reps of each rate in the log
## also collect the inter-gene collision UMI sequences for plotting
## NOTE: the umi_6 results need to be adj as we only have 4096 possibilities
## NOTE: need to divide the number of unique UMIs by the total number of unique
## UMIs allocated to find the obs/exp for vasa later

## theoretical
gather_collision_rates <- function(path) {
  ## first gather the directories
  analysis_dirs <- list.dirs(path = path,
                             full.names = TRUE,
                             recursive = FALSE)
  ## now dip into each and get the results
  sim_results <- lapply(analysis_dirs, function(x) {
    message("Working on dir: ", x)
    log_files <- list.files(x,
                            pattern = ".log",
                            full.names = TRUE)
    collision_rates <- lapply(log_files, function(y) {
      # parse the log file
      cr <- read.delim(y, header = FALSE)
      cr <- cr[grep("Inter-gene collision rate", cr$V1),]
      cr <- as.numeric(trimws(strsplit(cr, ":")[[1]][4]))
      return(cr)
    })
    collistion_rates <- do.call(c, collision_rates)
  })
  return(sim_results)
}

## VASA sim
umi_6_res <- gather_collision_rates(getwd())
names(umi_6_res) <- gsub("./", "", list.dirs(recursive = FALSE))
umi_6_res.long <- do.call(rbind, umi_6_res)
colnames(umi_6_res.long) <- paste0("rep_", 1:ncol(umi_6_res.long))
saveRDS(umi_6_res.long, "umi_6_simulation_results.rds")

## STORM sim
umi_8_res <- gather_collision_rates(getwd())
names(umi_8_res) <- gsub("./", "", list.dirs(recursive = FALSE))
umi_8_res.long <- do.call(rbind, umi_8_res)
colnames(umi_8_res.long) <- paste0("rep_", 1:ncol(umi_8_res.long))
saveRDS(umi_8_res.long, "umi_8_simulation_results.rds")

## ss3xpress sim
umi_10_res <- gather_collision_rates(getwd())
names(umi_10_res) <- gsub("./", "", list.dirs(recursive = FALSE))
umi_10_res.long <- do.call(rbind, umi_10_res)
colnames(umi_10_res.long) <- paste0("rep_", 1:ncol(umi_10_res.long))
saveRDS(umi_10_res.long, "umi_10_simulation_results.rds")

## sstotal umi sim
umi_16_res <- gather_collision_rates(getwd())
names(umi_16_res) <- gsub("./", "", list.dirs(recursive = FALSE))
umi_16_res.long <- do.call(rbind, umi_16_res)
colnames(umi_16_res.long) <- paste0("rep_", 1:ncol(umi_16_res.long))
saveRDS(umi_16_res.long, "umi_16_simulation_results.rds")

## pull in and plot
library(ggplot2)

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

# combine
all_sim_results <- rbind(umi_6_res.long,
                         umi_8_res.long,
                         umi_10_res.long,
                         umi_16_res.long)

ggplot(all_sim_results, aes(x = Var1, y = value, group = umi)) +
  geom_point(aes(color = umi)) +
  geom_smooth(method = "loess",
              se = FALSE,
              linetype = 2,
              color = "black") +
  xlab("Genes") +
  ylab("Inter-gene Colllision Rate") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
  scale_color_manual(values = c(`NNNNNNNN`="#63197FFF",
                                `NNNNNN`="#CF5917FF",
                                `NNNNNNNNNNNNNNNN` = "#8FC2FD",
                                `NNNNNNNNWW` = "#8CCC98"),
                     name = "UMI Base Diversity") +
  theme_bw(12) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
