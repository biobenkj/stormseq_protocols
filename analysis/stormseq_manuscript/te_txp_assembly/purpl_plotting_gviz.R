## GViz of LTR-PURPL for STORM
# NOTE: THIS IS MY CUSTOM GVIZ INSTALL
library(Gviz)
library(rtracklayer)

# Ensembl annotations
options(ucscChromosomeNames=FALSE)

# region
chrom <- "5"
from <- 27196384
to <- 27505749

# make a little granges
roi <- GRanges(seqnames = chrom,
               ranges = IRanges(start = from, end = to),
               strand = "*")

# ideogram track
ideoTrack <- IdeogramTrack(genome = "hg38",
                           chromosome = "chr5",
                           showTitle = FALSE,
                           showId = FALSE)
# not sure I really need this
# this is incomplete! can proceed if any primary chromosome
ideoTrack@bandTable$chrom <- gsub("chr|chrUn_",
                                  "", ideoTrack@bandTable$chrom)
ideoTrack@chromosome <- "5"

# genome axis track
gat <- GenomeAxisTrack(showTitle = FALSE,
                       labelPos = "above",
                       fontsize = 14)

# pull in the TE-txps and Purpl GTF
gtf <- import("purpl_gene_annots_with_tes_ens101.gtf")
purpl.gtf <- import("ens113.purpl.gtf")

# this is to aggregate all the short isoforms as there are many, many of them.
gtf.purpl <- purpl.gtf[purpl.gtf$gene_name %in% "PURPL",]
gtf.purpl.exon <- gtf.purpl[gtf.purpl$type %in% "exon",]
gtf.purpl.txp <- gtf.purpl[gtf.purpl$type %in% "transcript",]
table(width(gtf.purpl.txp) < 30000)
# FALSE  TRUE 
# 8   121

# pull out all the small transcripts and collapse
gtf.purpl.txp.short <- gtf.purpl.txp[width(gtf.purpl.txp) < 30000,]

# keep the transcripts
gtf.txp <- gtf[gtf$type %in% c("exon"),]
gtf.txp <- sort(gtf.txp)

# sep out the TE txps and PURPL txps
gtf.txp.purpl <- gtf.purpl.exon
gtf.txp.te <- gtf.txp[!gtf.txp$gene_name %in% "PURPL",]

# PURPL
gtf.txp.purpl.meta <- gtf.txp.purpl[gtf.txp.purpl$transcript_id %in% gtf.purpl.txp.short$transcript_id,]
gtf.txp.purpl <- gtf.txp.purpl[!gtf.txp.purpl$transcript_id %in% gtf.purpl.txp.short$transcript_id,]
gtf.txp.purpl.meta <- reduce(gtf.txp.purpl.meta)
# add back the needed mcols
mcols(gtf.txp.purpl.meta) <- data.frame(source = "havana_meta",
                                        type = "exon",
                                        score = as.numeric(NA),
                                        phase = as.integer(NA),
                                        gene_id = "ENSG00000250337",
                                        gene_name = "PURPL",
                                        transcript_id = "ENST00000250337",
                                        exon_number = c(1:length(gtf.txp.purpl.meta)),
                                        gene_version = "0",
                                        gene_source = "havana_meta",
                                        gene_biotype = "lncRNA",
                                        transcript_version = "0",
                                        transcript_name = "PURPL-000",
                                        transcript_source = "havana_meta",
                                        transcript_biotype = "lncRNA",
                                        tag = "basic",
                                        exon_id = paste0("ENSE000099999", 1:length(gtf.txp.purpl.meta)),
                                        exon_version = "1",
                                        transcript_support_level = as.character(NA))
gtf.txp.purpl <- c(gtf.txp.purpl,
                   gtf.txp.purpl.meta)

gtf.txp.purpl$transcript <- factor(gtf.txp.purpl$transcript_id)
gtf.txp.purpl$symbol <- factor(gtf.txp.purpl$gene_name)



# Now, create the GeneRegionTrack using the combined GRanges object.
purpl_annot_track <- GeneRegionTrack(gtf.txp.purpl,
                                     stacking = "squish",
                                     name = "TranscriptAnnotation",
                                     transcriptAnnotation = "symbol",
                                     showId = TRUE,
                                     geneSymbol = FALSE,
                                     showTitle = FALSE,
                                     fill = "#8405A7FF",
                                     col = "#8405A7FF",
                                     lwd = 2)

# TE-PURPL
gtf.txp.te$transcript <- factor(gtf.txp.te$transcript_id)
gtf.txp.te$symbol <- factor(ifelse(gtf.txp.te$transcript_id %in% "TU1199",
                                   "LTR1A2-PURPL.1",
                                   "LTR1A2-PURPL.2"))

te_purpl_annot_track <- GeneRegionTrack(gtf.txp.te,
                                        stacking = "squish",
                                        name = "TranscriptAnnotation",
                                        transcriptAnnotation = "symbol",
                                        showId = TRUE,
                                        geneSymbol = FALSE,
                                        showTitle = FALSE,
                                        fill = "black",
                                        col = "black",
                                        lwd = 2)

# check
# plotTracks(list(te_purpl_annot_track,
#                 purpl_annot_track))


# bring in the bams
# set the constant scale
cov_scale <- c(0,50)

# set the max value from the merged data
library(GenomicAlignments)
merged_bam <- readGAlignmentPairs(file = "merged_purpl_to_plot.bam",
                                  strandMode = 2)
merged_bam.junc <- summarizeJunctions(merged_bam)
max(merged_bam.junc$score)
# 227

rescale_scores <- function(scores) {
  fixed_max = 250
  (scores/fixed_max) * 100
}

# Purpl-unique txps
purpl_track <- AlignmentsTrack("purpl_known.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage","sashimi"),
                               sashimiScore = 2,
                               sashimiTransformation = rescale_scores)
displayPars(purpl_track) <- list(ylim = cov_scale,
                                 fill = "#8405A7FF",
                                 col.coverage = "#8405A7FF",
                                 col.sashimi = "#8405A7FF",
                                 showTitle = FALSE,
                                 stackHeight = 0.25)

# te-derived txp track
tetxp_track <- AlignmentsTrack("te_txp_only.bam",
                               isPaired = TRUE,
                               type = c("coverage","sashimi"),
                               sashimiScore = 2,
                               sashimiTransformation = rescale_scores)
displayPars(tetxp_track) <- list(ylim = cov_scale,
                                 fill = "#D35171FF",
                                 col.coverage = "#D35171FF",
                                 col.sashimi = "#D35171FF",
                                 showTitle = FALSE,
                                 stackHeight = 0.25)

# merged
merge_track <- AlignmentsTrack("merged_purpl_to_plot.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage", "sashimi"),
                               sashimiScore = 2,
                               sashimiTransformation = rescale_scores)

displayPars(merge_track) <- list(ylim = cov_scale,
                                 fill = "black",
                                 col.coverage = "black",
                                 col.sashimi = "black",
                                 showTitle = FALSE,
                                 stackHeight = 0.25)

# single cell with ltr1a2-purpl
sc_prp_track <- AlignmentsTrack("K3_purpl_te.bam",
                               isPaired = TRUE,
                               stacking = "squish",
                               type = c("coverage", "sashimi"),
                               sashimiScore = 1,
                               sashimiTransformation = rescale_scores)

displayPars(sc_prp_track) <- list(ylim = cov_scale,
                                 fill = "#FCA636FF",
                                 col.coverage = "#FCA636FF",
                                 col.sashimi = "#FCA636FF",
                                 showTitle = FALSE,
                                 stackHeight = 0.25)

# single cell withOUT ltr1a2-purpl
sc_nope_track <- AlignmentsTrack("E11_purpl_te.bam",
                                isPaired = TRUE,
                                stacking = "squish",
                                type = c("coverage", "sashimi"),
                                sashimiScore = 1,
                                sashimiTransformation = rescale_scores)

displayPars(sc_nope_track) <- list(ylim = cov_scale,
                                  fill = "#FECE91FF",
                                  col.coverage = "#FECE91FF",
                                  col.sashimi = "#FECE91FF",
                                  showTitle = FALSE,
                                  stackHeight = 0.25)


# bring in the TE annotations
te_annots <- import("repeatmasker_sorted.pints_filt.cut.bed")

teTrack <- AnnotationTrack(te_annots, name = "TE-track",
                           shape = "box", from = from, to = to,
                           fill = "black", col.line = "black",
                           showTitle = FALSE)

# build up the tracks to plot
plotTracks(list(ideoTrack, gat,
                merge_track, purpl_track, tetxp_track,
                sc_prp_track, sc_nope_track,
                teTrack,
                te_purpl_annot_track, purpl_annot_track),
           # showTitle = FALSE,
           chromosome = "5",
           from = from,
           to = to,
           #extend.left = 6000,
           col.axis = "black",
           col.title = "black", background.title = "white",
           cex.axis = 1, font.size = 14,
           sizes = c(0.2, 0.3,
                     0.7, 0.7, 0.7,
                     0.7, 0.7,
                     0.4,
                     0.2, 0.6))
