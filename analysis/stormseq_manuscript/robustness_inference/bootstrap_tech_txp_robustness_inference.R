## plot uncertainty estimates for HEK

## custom function using sleuth to import the per-cell
## abundance_1.h5 files

library(sleuth) # version 0.30.1

read_kb_bootstraps <- function(path = NULL,
                               file = "abundance_1.h5",
                               filter.nan = FALSE) {
  if (is.null(path)) {
    stop("Please provide a path to directory containing abundance files.")
  }
  
  if (getExtension(path) == "h5") {
    path <- gsub(paste0("/", file),
                 "", path)
  }
  
  file_to_read <- paste0(path, "/", file)
  if (!file.exists(file_to_read)) {
    stop("The file: ", file_to_read, " doesn't exist.")
  }
  
  # pre-flight checks complete
  # read in file
  message("Reading in: ", file_to_read)
  kb <- sleuth:::read_kallisto(file_to_read,
                               read_bootstrap = TRUE)
  # convert bootstraps to a matrix
  message("Summarizing bootstraps...")
  kb.summary <- suppressWarnings(sleuth:::summarize_bootstrap(kb))
  kb.summary.df <- as.data.frame(kb.summary)
  
  # filter
  if (filter.nan) {
    message("Filtering out NaN entries for bs_cv_tpm...")
    kb.summary.df <- kb.summary.df[!is.nan(kb.summary.df$bs_cv_tpm),]
  }
  message("Done.")
  return(kb.summary.df)
}

# helper
getExtension <- function(file){ 
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
} 

summarize_within_subsamples <- function(df) {
  # Split the data by target_id
  split_data <- split(df, df$target_id)
  
  # Function to summarize metrics for a single target_id
  summarize_single_target_id <- function(df_subset) {
    # Calculate detection rate: proportion of subsamples where the transcript is detected
    detection_rate <- sum(df_subset$bs_mean_tpm > 0) / nrow(df_subset)
    
    # Weighted mean expression considering detection rate
    weighted_mean_expression <- mean(df_subset$bs_mean_tpm * detection_rate, na.rm = TRUE)
    
    # Calculate adjusted CV using detection rate
    adjusted_cv <- sd(df_subset$bs_mean_tpm, na.rm = TRUE) / weighted_mean_expression
    
    # Summarize metrics
    summarized_df <- data.frame(
      target_id = unique(df_subset$target_id),
      mean_of_means = mean(df_subset$bs_mean_tpm, na.rm = TRUE),
      sd_of_means = sd(df_subset$bs_mean_tpm, na.rm = TRUE),
      var_of_means = var(df_subset$bs_mean_tpm, na.rm = TRUE),
      mean_of_sds = mean(df_subset$bs_sd_tpm, na.rm = TRUE),
      mean_cv_across_subsamples = mean(df_subset$bs_cv_tpm, na.rm = TRUE),
      median_cv_across_subsamples = median(df_subset$bs_cv_tpm, na.rm = TRUE),
      detection_rate = detection_rate,
      weighted_mean_expression = weighted_mean_expression,
      adjusted_cv = adjusted_cv
    )
    
    return(summarized_df)
  }
  
  # Apply the summarization function to each target_id
  summarized_list <- lapply(split_data, summarize_single_target_id)
  
  # Collate the results back together into a single data frame
  summarized_df <- do.call(rbind, summarized_list)
  
  return(summarized_df)
}

calculate_subsample_metrics <- function(summarized_df_list) {
  # Function to calculate pooled standard deviation and other metrics across subsamples
  
  # Extract the required columns and stack them together
  stacked_data <- do.call(rbind, summarized_df_list)
  
  # Split the combined data by target_id
  split_data <- split(stacked_data, stacked_data$target_id)
  
  # Group by target_id to calculate metrics across subsamples
  summarize_across_subsamples <- function(df_subset) {
    # Convert relevant columns to numeric, in case they were converted to characters
    df_subset$mean_of_means <- as.numeric(df_subset$mean_of_means)
    df_subset$sd_of_means <- as.numeric(df_subset$sd_of_means)
    df_subset$var_of_means <- as.numeric(df_subset$var_of_means)
    df_subset$mean_of_sds <- as.numeric(df_subset$mean_of_sds)
    df_subset$mean_cv_across_subsamples <- as.numeric(df_subset$mean_cv_across_subsamples)
    df_subset$median_cv_across_subsamples <- as.numeric(df_subset$median_cv_across_subsamples)
    df_subset$detection_rate <- as.numeric(df_subset$detection_rate)
    df_subset$weighted_mean_expression <- as.numeric(df_subset$weighted_mean_expression)
    df_subset$adjusted_cv <- as.numeric(df_subset$adjusted_cv)
    
    # Summarize metrics
    summarized_df <- data.frame(
      target_id = unique(df_subset$target_id),
      mean_of_means = mean(df_subset$mean_of_means, na.rm = TRUE),
      sd_of_means = mean(df_subset$sd_of_means, na.rm = TRUE),  # Mean of sd_of_means across subsamples
      var_of_means = mean(df_subset$var_of_means, na.rm = TRUE),  # Mean of var_of_means across subsamples
      mean_of_sds = mean(df_subset$mean_of_sds, na.rm = TRUE),
      mean_cv_across_subsamples = mean(df_subset$mean_cv_across_subsamples, na.rm = TRUE),
      median_cv_across_subsamples = median(df_subset$median_cv_across_subsamples, na.rm = TRUE),
      mean_detection_rate = mean(df_subset$detection_rate, na.rm = TRUE),
      median_detection_rate = median(df_subset$detection_rate, na.rm = TRUE),
      mean_weighted_expression = mean(df_subset$weighted_mean_expression, na.rm = TRUE),
      sd_weighted_expression = sd(df_subset$weighted_mean_expression, na.rm = TRUE),
      mean_adjusted_cv = mean(df_subset$adjusted_cv, na.rm = TRUE),
      median_adjusted_cv = median(df_subset$adjusted_cv, na.rm = TRUE)
    )
    
    return(summarized_df)
  }
  
  # Apply the summarization function to each target_id
  summarized_list <- lapply(split_data, summarize_across_subsamples)
  
  # Collate the results back together into a single data frame
  final_metrics_summary <- do.call(rbind, summarized_list)
  
  return(final_metrics_summary)
}

# need this for mclapply()
library(parallel)

## STORM
# sample01-sample10
# do each per random seed
samp_dirs <- "analysis/bootstrap/storm/sample10"
abund_dirs <- list.files(samp_dirs,
                         pattern = "abundance_1.h5",
                         full.names = TRUE,
                         recursive = TRUE)
# this is best done in a high performance computing environment
boot_quants.storm <- mclapply(abund_dirs, function(x) {
  return(read_kb_bootstraps(x))
}, mc.cores = 40)

# load libraries if starting here
library(dplyr)
library(tidyr)
# read in the sleuth processed .rds file per random seed
boot_quants.storm <- readRDS("storm_sample10_raw_bootstraps.rds")

# collate
hek_boot_quants.storm.collated <- do.call(rbind, boot_quants.storm)

# summarize
hek.storm.summarized <- summarize_within_subsamples(hek_boot_quants.storm.collated)

# save
saveRDS(hek.storm.summarized, file = "storm_sample10_raw_bootstraps_summarized.rds")



## VASA
# sample01-sample10
samp_dirs <- "analysis/bootstrap/vasa/sample10"
abund_dirs <- list.files(samp_dirs,
                         pattern = "abundance_1.h5",
                         full.names = TRUE,
                         recursive = TRUE)
# read in data on HPC
boot_quants.vasa <- mclapply(abund_dirs, function(x) {
  return(read_kb_bootstraps(x))
}, mc.cores = 40)

# read in the sleuth processed .rds file per random seed
boot_quants.vasa <- readRDS("vasa_sample10_raw_bootstraps.rds")

# collate
hek_boot_quants.vasa.collated <- do.call(rbind, boot_quants.vasa)

# summarize
hek.vasa.summarized <- summarize_within_subsamples(hek_boot_quants.vasa.collated)
saveRDS(hek.vasa.summarized, file = "vasa_sample10_raw_bootstraps_summarized.rds")


## SSTOTAL
# sample01-sample10
samp_dirs <- "analysis/bootstrap/sstotal/sample10"
abund_dirs <- list.files(samp_dirs,
                         pattern = "abundance_1.h5",
                         full.names = TRUE,
                         recursive = TRUE)
# read in data on HPC
boot_quants.sstotal <- mclapply(abund_dirs, function(x) {
  return(read_kb_bootstraps(x))
}, mc.cores = 40)

# read in the sleuth processed .rds file per random seed
boot_quants.sstotal <- readRDS("sstotal_sample10_raw_bootstraps.rds")

# collate
hek_boot_quants.sstotal.collated <- do.call(rbind, boot_quants.sstotal)

# summarize
hek.sstotal.summarized <- summarize_within_subsamples(hek_boot_quants.sstotal.collated)
saveRDS(hek.sstotal.summarized, file = "sstotal_sample10_raw_bootstraps_summarized.rds")

# directory names too to pull in per-tech summarized rds files
tech <- c("storm", "vasa",
          "sstotal")

# load libraries if starting here
library(dplyr)
library(tidyr)
all_tech_sds <- lapply(tech, function(x) {
  # Read in the pooled sds
  tech_pooled_sd <- lapply(list.files("/path/to/summarized/rds_files/",
                                      pattern = x,
                                      full.names = TRUE), function(y) {
                                        pooled_in <- as.matrix(readRDS(y))
                                        # Convert to data frame and ensure the row names are a column for joining
                                        return(as.data.frame(pooled_in))
                                      })
  
  # calculate the subsample metrics using the above functions
  tech_pooled_summary <- calculate_subsample_metrics(tech_pooled_sd)
  
  return(tech_pooled_summary)
})

names(all_tech_sds) <- tech

# remove any NAs from each
all_tech_sds <- lapply(all_tech_sds, function(x) {
  return(x[complete.cases(x),])
})

# filter to transcripts that have at least 10% mean detection rates across subsamples
all_tech_sds.filt <- lapply(all_tech_sds, function(x) {
  return(x[x$mean_detection_rate >= 0.1,])
})

# find the intersection of transcripts across all technologies
extract_target_ids <- function(df) {
  return(df$target_id)
}
target_ids_list <- lapply(all_tech_sds.filt, extract_target_ids)

# find cross technology common txps
common_transcripts <- Reduce(intersect, target_ids_list)

all_tech_sds.filt <- lapply(all_tech_sds.filt, function(x) {
  return(x[x$target_id %in% common_transcripts,])
})

tech_summaries <- all_tech_sds.filt

tech_names <- c("STORM-seq", "VASA-seq",
                "Smart-seq-total")
for (i in seq_along(tech_summaries)) {
  tech_summaries[[i]]$technology <- tech_names[i]
}

combined_summary <- bind_rows(tech_summaries)

combined_summary$technology <- factor(combined_summary$technology,
                                      levels = c("Smart-seq-total", "VASA-seq",
                                                 "STORM-seq"))

# import the Ensembl 101 GTF
gtf <- rtracklayer::import("Homo_sapiens.GRCh38.101.gtf")

# Extract transcript lengths
transcript_lengths <- gtf %>%
  as.data.frame() %>%
  filter(type == "transcript") %>%
  mutate(length = end - start) %>%
  select(transcript_id = transcript_id, length)

# Define bins
breaks <- c(0, 500, 2000, 4000, Inf)
labels <- c("0-500", "500-2000", "2000-4000", ">4000")

# Create bins
exon_lengths <- gtf %>%
  as.data.frame() %>%
  filter(type == "exon") %>%
  mutate(length = end - start + 1) %>%
  group_by(transcript_id) %>%
  summarise(cDNA_length = sum(length))

exon_lengths <- exon_lengths %>%
  mutate(length_bin = cut(cDNA_length, breaks = breaks, labels = labels, include.lowest = TRUE))

transcript_lengths <- transcript_lengths %>%
  mutate(length_bin = cut(length, breaks = breaks, labels = labels, include.lowest = TRUE))

# fix the gene rownames
fix_rownames <- function(txis, stub="^ENS", sep="\\.", idx=1) { 
  
  fixable <- grep(stub, txis$target_id)
  if (length(fixable) < 1) { 
    message("No fixable rownames. Returning unaltered.")
    return(txis)
  }
  txis$target_id[fixable] <- 
    sapply(strsplit(txis$target_id[fixable], sep), `[`, idx)
  message("Fixed ", length(fixable), " rownames.") 
  return(txis) 
  
}
combined_summary <- fix_rownames(combined_summary)

# Merge transcript lengths with calculated data
combined_summary <- combined_summary %>%
  left_join(exon_lengths, by = c("target_id" = "transcript_id"))

# remove the 0-500 bin since they are often cleaned away during library prep...
combined_summary.filt <- combined_summary[!combined_summary$length_bin %in% "0-500",]
combined_summary.filt <- combined_summary.filt[complete.cases(combined_summary.filt),]
combined_summary.filt$length_bin <- factor(combined_summary.filt$length_bin,
                                           levels = c("500-2000", "2000-4000",
                                                      ">4000"))
library(ggplot2)
library(viridis)

# Plotting detection rate values across technologies for each bin
ggplot(combined_summary.filt, aes(y = mean_detection_rate, x = technology, fill = technology)) +
  geom_violin() +
  stat_summary(fun = median, geom = "crossbar", size = 0.8) +
  facet_wrap(~ length_bin, scales = "free_y") +
  scale_fill_manual(values = c(`STORM-seq`="#63197FFF",
                               `VASA-seq`="#CF5917FF",
                               `Smart-seq-total` = "#8FC2FD"),
                    name = "Technology") +
  ggtitle("Transcript Detection Rate (Higher is Better)") +
  ylab("Bootstrapped Detection Rate") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", linewidth = 0.8) +
  theme_bw(12) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# variance/robustness plot
ggplot(combined_summary.filt, aes(x = technology, y = log2(mean_adjusted_cv + 1), fill = technology)) +
  geom_violin() +
  stat_summary(fun = median, geom = "crossbar", size = 0.8) +
  scale_fill_manual(values = c(`STORM-seq`="#63197FFF",
                               `VASA-seq`="#CF5917FF",
                               `Smart-seq-total` = "#8FC2FD"),
                    name = "Technology") +
  facet_wrap(~ length_bin, scales = "free_y") +
  ggtitle("Transcript Detection Robustness (Lower is Better)") +
  ylab("Bootstrapped Transcript Variance") +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8,2)) +
  theme_bw(12) +
  theme(
    plot.title = element_text(size = 18,
                              hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
