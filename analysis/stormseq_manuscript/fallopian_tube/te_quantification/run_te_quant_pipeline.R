## code to run the FTE TE quantification pipeline
## source annotation was Ensembl 101
## repeat masker annotation hg38 - Dec 2013 - RepeatMasker open-4.0.5 - Repeat Library 20140131


## source the TE quantification pipeline R script
source("/varidata/research/projects/shen/projects/2024_Bens_swan_song/src/TE_quantification_pipeline/TE_quantification_pipeline.R")
# data.table 1.14.8 using 40 threads (see ?getDTthreads).  Latest news: r-datatable.com
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:data.table’:
#   
#   between, first, last
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union

# set up parameters
# NOTE: this is for patient 1 FTE
target_dir <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/fte/pat1/te_quant"
config_file <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/src/TE_quantification_pipeline/supporting_files/config_file.tsv"
read_names_file <- "/varidata/research/projects/shen/projects/2024_Bens_swan_song/src/storm_fte/te_quant/pat1/read_names.tsv"
sample <- "STORM_FTE_pat1"
create_count_matrix <- "TRUE"
scuttle_filter <- "TRUE"
log_enrich <- "TRUE"

# run it
quantify_TEs(config_file,
             read_names_file,
             sample,
             target_dir,
             create_count_matrix,
             scuttle_filter,
             log_enrich)

## session info
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-conda-linux-gnu (64-bit)
# Running under: AlmaLinux 9.5 (Teal Serval)
# 
# Matrix products: default
# BLAS/LAPACK: /varidata/research/projects/shen/tools/benkj/mambaforge/envs/te/lib/libopenblasp-r0.3.23.so
# 
# locale:
#   [1] LC_CTYPE=C.utf8       LC_NUMERIC=C          LC_TIME=C.utf8
# [4] LC_COLLATE=C.utf8     LC_MONETARY=C.utf8    LC_MESSAGES=C.utf8
# [7] LC_PAPER=C.utf8       LC_NAME=C             LC_ADDRESS=C
# [10] LC_TELEPHONE=C        LC_MEASUREMENT=C.utf8 LC_IDENTIFICATION=C
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] stringr_1.5.0     scales_1.2.1      ggplot2_3.4.3     dplyr_1.1.3
# [5] data.table_1.14.8
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.4      withr_2.5.1      utf8_1.2.3       grid_4.2.3
# [5] R6_2.5.1         gtable_0.3.4     lifecycle_1.0.3  magrittr_2.0.3
# [9] pillar_1.9.0     stringi_1.7.12   rlang_1.1.1      cli_3.6.1
# [13] vctrs_0.6.3      generics_0.1.3   tools_4.2.3      glue_1.6.2
# [17] munsell_0.5.0    compiler_4.2.3   colorspace_2.1-0 pkgconfig_2.0.3
# [21] tidyselect_1.2.0 tibble_3.2.1
