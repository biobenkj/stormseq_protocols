# plot the eRNA screenshot Hui wants

library(Gviz)
library(rtracklayer)

# Ensembl annotations
options(ucscChromosomeNames=FALSE)

# add some white space
# options(Gviz.leftMargin = 300)

# region
chrom <- "11"
from <- 61971495
to <- 61974026

# make a little granges
roi <- GRanges(seqnames = chrom,
               ranges = IRanges(start = from, end = to),
               strand = "*")

# ideogram track
ideoTrack <- IdeogramTrack(genome = "hg38",
                           chromosome = "chr11",
                           showTitle = FALSE,
                           showId = FALSE)
# not sure I really need this
# this is incomplete! can proceed if any primary chromosome
ideoTrack@bandTable$chrom <- gsub("chr|chrUn_",
                                  "", ideoTrack@bandTable$chrom)
ideoTrack@chromosome <- "11"

# genome axis track
gat <- GenomeAxisTrack(showTitle = FALSE,
                       labelPos = "above",
                       fontsize = 14)

# add the enhancers
enh_track <- rtracklayer::import("~/Documents/manuscripts/storm_seq/erna/Distal_K562_gro_pro_cap_pints_hg38_gencode24.only_regions.bed")
seqlevelsStyle(enh_track) <- "Ensembl"

enhTrack <- AnnotationTrack(enh_track, name = "Enhancers",
                              shape = "box", from = from, to = to,
                              fill = "gray", showTitle = FALSE)

# bring in the bams
# set the constant scale
cov_scale <- c(0,100)

# cell 1
cell1_track <- AlignmentsTrack("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/plots/igv_shot_chr11/I17_soloAligned.sortedByCoord.out.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage","pileup"))
displayPars(cell1_track) <- list(ylim = cov_scale,
                                 fill = "#8405A7FF",
                                 col.coverage = "#8405A7FF",
                                 fill.reads = "#8405A7FF",
                                 showTitle = FALSE,
                                 stackHeight = 0.5)

cell2_track <- AlignmentsTrack("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/plots/igv_shot_chr11/D14_soloAligned.sortedByCoord.out.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage","pileup"))
displayPars(cell2_track) <- list(ylim = cov_scale,
                                 fill = "#D35171FF",
                                 col.coverage = "#D35171FF",
                                 fill.reads = "#D35171FF",
                                 showTitle = FALSE,
                                 stackHeight = 0.5)

cell3_track <- AlignmentsTrack("~/Documents/manuscripts/storm_seq/erna/k562_storm/storm/plots/igv_shot_chr11/J5_soloAligned.sortedByCoord.out.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage","pileup"))

displayPars(cell3_track) <- list(ylim = cov_scale,
                                 fill = "black",
                                 col.coverage = "black",
                                 fill.reads = "black",
                                 showTitle = FALSE,
                                 stackHeight = 0.5)

# bring in the tt-seq fwd and rev bigwigs for back up

rep1_ttseq <- rtracklayer::import("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/msb_2021/GSM4610686_L_K562_Rep1.coverage.track.combined_fwd_rev.no_chr.bw")
ttseqtrack <- DataTrack(rep1_ttseq,
                       type = c("polygon"), chromosome = "11",
                       from = from, to = to,
                       fill.mountain = c("gray", "#FECE91FF"),
                       col.mountain = "transparent", name = "TT-seq",
                       showTitle = F, stackHeight = 0.5,
                       ylim = cov_scale)
rep2_ttseq <- rtracklayer::import("~/Documents/manuscripts/storm_seq/erna/k562_storm/tt_seq/msb_2021/GSM4610687_L_K562_Rep2.coverage.track.combined_fwd_rev.no_chr.bw")
ttseqtrack2 <- DataTrack(rep2_ttseq,
                        type = c("polygon"), chromosome = "11",
                        from = from, to = to,
                        fill.mountain = c("gray", "#2A788EFF"),
                        col.mountain = "transparent", name = "TT-seq",
                        showTitle = F, stackHeight = 0.5,
                        ylim = cov_scale)


# build up the tracks to plot
plotTracks(list(ideoTrack, gat,
                enhTrack,
                cell1_track, cell2_track, cell3_track,
                ttseqtrack, ttseqtrack2),
           # showTitle = FALSE,
           chromosome = "11",
           from = from,
           to = to,
           #extend.left = 6000,
           col.axis = "black",
           col.title = "black", background.title = "white",
           cex.axis = 1, font.size = 14,
           sizes = c(0.2, 0.3,
                     0.1,
                     1, 1, 0.6,
                     0.25, 0.25))

