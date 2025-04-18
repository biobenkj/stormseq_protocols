## look at where the mismapped reads/outside annotations are going

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome)

# R-loop forming sequence retrieval from RLHub
library(RLHub)
library(ExperimentHub)
eh <- ExperimentHub()
resource_id <- "EH6806" # this is the metadata which is the only one that contains ranges
rl_hg38 <- eh[[resource_id]]
# extract regions that overlap with R loop forming sequences
rl_hg38_is_rlfs <- rl_hg38[rl_hg38$is_rlfs,]
# convert to granges
rl_hg38_is_rlfs.gr <- as(gsub(":\\.", "", rl_hg38_is_rlfs$location),
                         "GRanges")
seqlevelsStyle(rl_hg38_is_rlfs.gr) <- "Ensembl"

rl_hg38.gr <- as(gsub(":\\.", "", rl_hg38$location),
                 "GRanges")
seqlevelsStyle(rl_hg38.gr) <- "Ensembl"

# actual genome used for alignment
genome <- readRDS("Homo_sapiens.GRCh38.dna.primary_assembly_ercc92_dnastringset.rds")

# exclude ERCCs
genome <- genome[grep("ERCC-", invert = TRUE, names(genome))]

# import the polyA and polyT runs bed files
polyA_runs <- import("hg38_polyA_min6_c2_at_filt.bed.gz")
polyT_runs <- import("hg38_polyT_min6_c2_at_filt.bed.gz")

# Function to deduplicate ranges based on cell barcode
library(stringr)
process_granges <- function(gr) {
  # Ensure the required package is loaded
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("The 'stringr' package is required but not installed. Please install it using install.packages('stringr').")
  }
  
  # Check if names are present
  if (is.null(gr$name)) {
    stop("The GRanges object must have names to extract the 'CB:' field.")
  }
  
  # Extract names
  gr_names <- gr$name
  
  # Check for 'CB:' entries
  has_cb <- str_detect(gr_names, "CB:")
  
  # Initialize result vector
  cb_values <- rep(NA_character_, length(gr_names))
  
  # Extract 'CB:' substrings
  cb_substrings <- str_extract(gr_names[has_cb], "CB:[^;]+")
  
  # Remove 'CB:' prefix to get the barcode
  cb_values[has_cb] <- str_remove(cb_substrings, "CB:")
  
  # Add 'CB' as a metadata column
  mcols(gr)$CB <- cb_values
  
  # Create a data.frame with relevant columns for duplication checking
  df <- data.frame(
    seqnames = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr),
    strand = as.character(strand(gr)),
    CB = mcols(gr)$CB,
    stringsAsFactors = FALSE
  )
  
  # Identify duplicates based on CB and exact ranges
  dup_indices <- duplicated(df)
  
  # Subset the GRanges object to remove duplicates
  gr_unique <- gr[!dup_indices]
  
  # Return the de-duplicated GRanges object
  return(gr_unique)
}

# function to calculate distance to a polyN run
compute_polyA_polyT_distances <- function(my_gr, polyA_ranges, polyT_ranges, threshold) {
  # Ensure that the required package is loaded
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package is required but not installed. Please install it using install.packages('GenomicRanges').")
  }
  
  # Check that all inputs are GRanges objects
  if (!inherits(my_gr, "GRanges") || !inherits(polyA_ranges, "GRanges") || !inherits(polyT_ranges, "GRanges")) {
    stop("All inputs must be GRanges objects.")
  }
  
  # Compute distance to nearest polyA ranges
  nearest_polyA <- distanceToNearest(my_gr, polyA_ranges, ignore.strand = TRUE)
  
  # Initialize distance vector with NA
  polyA_distances <- rep(NA_integer_, length(my_gr))
  
  # Assign distances
  polyA_distances[queryHits(nearest_polyA)] <- mcols(nearest_polyA)$distance
  
  # Compute distance to nearest polyT ranges
  nearest_polyT <- distanceToNearest(my_gr, polyT_ranges, ignore.strand = TRUE)
  
  # Initialize distance vector with NA
  polyT_distances <- rep(NA_integer_, length(my_gr))
  
  # Assign distances
  polyT_distances[queryHits(nearest_polyT)] <- mcols(nearest_polyT)$distance
  
  # Add distances to mcols
  mcols(my_gr)$polyA_dist <- polyA_distances
  mcols(my_gr)$polyT_dist <- polyT_distances
  
  # Add new logical column based on threshold
  # Check if either polyA_dist or polyT_dist is less than the threshold
  # Handle NA values appropriately
  mcols(my_gr)$within_threshold <- with(mcols(my_gr), 
                                        (polyA_dist <= threshold | polyT_dist <= threshold))
  
  # Replace NA in 'within_threshold' with FALSE
  mcols(my_gr)$within_threshold[is.na(mcols(my_gr)$within_threshold)] <- FALSE
  
  # Return the updated GRanges object
  return(my_gr)
}

count_barcodes <- function(strings, has_tag = TRUE) {
  if (has_tag) {
    # Process strings with 'CB:' tags
    parts_list <- strsplit(strings, ";")
    
    cb_values <- sapply(parts_list, function(parts) {
      cb_part <- parts[grepl("^CB:", parts)]
      if (length(cb_part) > 0) {
        sub("^CB:", "", cb_part)
      } else {
        NA
      }
    })
    
    # Remove any NA values
    cb_values <- cb_values[!is.na(cb_values)]
    
  } else {
    # Process strings without 'CB:' tags; extract up to the first '.'
    cb_values <- sapply(strings, function(s) {
      dot_pos <- regexpr("\\.", s)
      if (dot_pos > 0) {
        substr(s, 1, dot_pos - 1)
      } else {
        # If no '.', use the whole string
        s
      }
    })
  }
  
  # Count the occurrences of each barcode
  barcode_counts <- table(cb_values)
  
  # Convert the table to a data frame
  barcode_counts_df <- as.data.frame(barcode_counts)
  colnames(barcode_counts_df) <- c("Cell_Barcode", "Count")
  
  return(barcode_counts_df)
}

## vasa
vasa_reads <- rtracklayer::import("merged_hek293t_nonannot_reads.sorted.bed")

# exclude ERCCs
vasa_reads <- vasa_reads[grep("ERCC-", invert = TRUE, seqnames(vasa_reads)),]
seqlevels(vasa_reads) <- seqlevels(vasa_reads)[grep("ERCC-", invert = TRUE, seqlevels(vasa_reads))]

# add seqlengths
seqlengths(vasa_reads) <- seqlengths(genome)[names(seqlengths(vasa_reads))]

# count mispriming within cells,
# filtering for R-loop regions
vasa_reads.filt.rloopfilt <- subsetByOverlaps(vasa_reads,
                                              rl_hg38.gr,
                                              invert = TRUE)
vasa_mispriming_events_percell <- count_barcodes(vasa_reads.filt.rloopfilt$name)


## sstotal
sstotal_reads <- rtracklayer::import("~/Documents/manuscripts/storm_seq/genomic_region_analysis/sstotal_umi/merged_hek293t_sstotal_mismapped.sorted.bed")

# exclude ERCCs
sstotal_reads <- sstotal_reads[grep("ERCC-", invert = TRUE, seqnames(sstotal_reads)),]
seqlevels(sstotal_reads) <- seqlevels(sstotal_reads)[grep("ERCC-", invert = TRUE, seqlevels(sstotal_reads))]

# add seqlengths
seqlengths(sstotal_reads) <- seqlengths(genome)[names(seqlengths(sstotal_reads))]

# count mispriming within cells,
# filtering for R-loop regions
sstotal_reads.filt.rloopfilt <- subsetByOverlaps(sstotal_reads,
                                              rl_hg38.gr,
                                              invert = TRUE)
sstotal_mispriming_events_percell <- count_barcodes(sstotal_reads.filt.rloopfilt$name,
                                                    has_tag = FALSE)


# bulk
bulk_reads <- import("~/Documents/manuscripts/storm_seq/genomic_region_analysis/bulk_total/merged_bulk_100k_ds_tworep_nonannot.sorted.converted_frags.bed")

# exclude ERCCs
bulk_reads <- bulk_reads[grep("ERCC-", invert = TRUE, seqnames(bulk_reads)),]
seqlevels(bulk_reads) <- seqlevels(bulk_reads)[grep("ERCC-", invert = TRUE, seqlevels(bulk_reads))]

# add seqlengths
seqlengths(bulk_reads) <- seqlengths(genome)[names(seqlengths(bulk_reads))]

# count mispriming within cells,
# filtering for R-loop regions
bulk_reads.filt.rloopfilt <- subsetByOverlaps(bulk_reads,
                                                 rl_hg38.gr,
                                                 invert = TRUE)

bulk_mispriming_events_percell <- count_barcodes(bulk_reads.filt.rloopfilt$name,
                                                 has_tag = FALSE)

## ss3xpress
ss3_reads <- rtracklayer::import("~/Documents/manuscripts/storm_seq/genomic_region_analysis/ss3/merged_nonannot_reads_ss3xpress_ens101.sorted.converted_frags.bed")

# exclude ERCCs
ss3_reads <- ss3_reads[grep("ERCC-", invert = TRUE, seqnames(ss3_reads)),]
seqlevels(ss3_reads) <- seqlevels(ss3_reads)[grep("ERCC-", invert = TRUE, seqlevels(ss3_reads))]

# add seqlengths
seqlengths(ss3_reads) <- seqlengths(genome)[names(seqlengths(ss3_reads))]

# fix the cell barcode names
ss3_reads$name <- gsub("Aligned", "", ss3_reads$name)

# count mispriming within cells,
# filtering for R-loop regions
ss3_reads.filt.rloopfilt <- subsetByOverlaps(ss3_reads,
                                             rl_hg38.gr,
                                             invert = TRUE)

ss3_mispriming_events_percell <- count_barcodes(ss3_reads$name,
                                                has_tag = FALSE)


# storm
storm_reads <- rtracklayer::import("~/Documents/manuscripts/storm_seq/genomic_region_analysis/storm/merged_hek293t_storm_nonannot_reads_ens101.sorted.converted_frags.bed")

# exclude ERCCs
storm_reads <- storm_reads[grep("ERCC-", invert = TRUE, seqnames(storm_reads)),]
seqlevels(storm_reads) <- seqlevels(storm_reads)[grep("ERCC-", invert = TRUE, seqlevels(storm_reads))]

# add seqlengths
seqlengths(storm_reads) <- seqlengths(genome)[names(seqlengths(storm_reads))]

# fix the cell barcode names
storm_reads$name <- gsub("_soloAligned", "", storm_reads$name)

# count mispriming within cells,
# filtering for R-loop regions
storm_reads.filt.rloopfilt <- subsetByOverlaps(storm_reads,
                                               rl_hg38.gr,
                                               invert = TRUE)

storm_mispriming_events_percell <- count_barcodes(storm_reads$name,
                                                has_tag = FALSE)


# plot a tie figher plot of mispriming events per cell
vasa_mispriming_events_percell$technology <- "VASA-seq"
sstotal_mispriming_events_percell$technology <- "Smart-seq-total"
storm_mispriming_events_percell$technology <- "STORM-seq"
ss3_mispriming_events_percell$technology <- "Smart-seq3xpress"
bulk_mispriming_events_percell$technology <- "Bulk"
all_mispriming <- rbind(bulk_mispriming_events_percell,
                        storm_mispriming_events_percell,
                        vasa_mispriming_events_percell,
                        sstotal_mispriming_events_percell,
                        ss3_mispriming_events_percell)

library(dplyr)

# Calculate mean and standard deviation by technology
summary_data <- all_mispriming %>%
  group_by(technology) %>%
  summarise(
    mean_count = mean(Count),
    sd_count = sd(Count)
  )

# Calculate the mean of "Bulk"
bulk_mean <- summary_data %>%
  filter(technology == "Bulk") %>%
  pull(mean_count)

# Subtract the "Bulk" mean from all means
summary_data <- summary_data %>%
  mutate(
    adjusted_mean = mean_count - bulk_mean
  )

summary_data$technology <- factor(summary_data$technology,
                                  levels = c("Bulk",
                                             "Smart-seq3xpress",
                                             "STORM-seq",
                                             "Smart-seq-total",
                                             "VASA-seq"))

# percentages
summary_data$percent_total <- (summary_data$mean_count/100000) * 100
summary_data$percent_total_sd <- (summary_data$sd_count/100000) * 100

library(ggplot2)

ggplot(summary_data, aes(x = technology, y = percent_total, color = technology)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = percent_total - percent_total_sd, ymax = percent_total + percent_total_sd),
                width = 0.4, size = 1) +
  ylab("Mean Percent ± SD") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_color_manual(values = c(`STORM-seq`="#63197FFF",
                                `VASA-seq`="#CF5917FF",
                                `Smart-seq-total` = "#8FC2FD",
                                `Smart-seq3xpress` = "#8CCC98",
                                `Bulk` = "goldenrod"),
                     name = "Technology") +
  ggtitle("Per-cell % reads aligning to genomic background",
          subtitle = "100k reads/cell") +
  coord_flip() +
  theme_bw(12) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    plot.title = element_text(hjust = 0.5,
                              size = 18),
    plot.subtitle = element_text(hjust = 0.5,
                                 size = 16))
