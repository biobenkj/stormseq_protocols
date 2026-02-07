## code to create the read_names.tsv file

# list the file prefixes
storm_bams <- list.files(path = "/varidata/research/projects/shen/projects/2024_Bens_swan_song/analysis/te_quant/storm/fte/pat2/te_quant/aligned_files",
                         pattern = "*Aligned.sortedByCoord.out.bam$")
storm_cells <- gsub("Aligned.sortedByCoord.out.bam",
                    "",
                    storm_bams)

# other pieces to skip pipeline steps
# and reconstruct the read_names.tsv file
rnames <- data.frame(cells = storm_cells,
                     paired = rep("TRUE", length(storm_cells)),
                     raw_reads_1 = rep("NULL", length(storm_cells)),
                     raw_reads_2 = rep("NULL", length(storm_cells)),
                     trim = rep("FALSE", length(storm_cells)),
                     trimmed_reads_1 = rep("NULL", length(storm_cells)),
                     trimmed_reads_2 = rep("NULL", length(storm_cells)),
                     align = rep("FALSE", length(storm_cells)),
                     aligned_files = storm_bams,
                     fc = rep("TRUE", length(storm_cells)),
                     fc_files = rep("NULL", length(storm_cells)),
                     sc = rep("TRUE", length(storm_cells)))

# write it out
write.table(rnames,
            file = "read_names.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')