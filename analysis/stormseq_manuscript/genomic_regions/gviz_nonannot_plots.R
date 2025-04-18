## now for the plotting
library(Gviz)

# region
chrom <- 6
from <- 28896187
to <- 28896388

# make a little granges
roi <- GRanges(seqnames = "chr6",
               ranges = IRanges(start = from, end = to),
               strand = "*")

# ideogram track
ideoTrack <- IdeogramTrack(genome = "hg38",
                           chromosome = chrom,
                           showTitle = FALSE,
                           showId = FALSE)

# genome axis track
gat <- GenomeAxisTrack(showTitle = FALSE,
                       labelPos = "below", fontsize = 18)

# bring in the gene annotations
library(biomaRt)
# Use the Ensembl 101 host
bm <- useEnsembl(host = "https://nov2020.archive.ensembl.org", 
                 biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "hg38", chromosome = chrom, 
                                    start = from, end = to,
                                    name = "ENSEMBL", biomart = bm,
                                    collapseTranscripts = "meta",
                                    stacking = "squish",
                                    showTitle = F, fill = "black")

# bring in the sequence context as boxes
genome <- readRDS("Homo_sapiens.GRCh38.dna.primary_assembly_ercc92_dnastringset.rds")
names(genome) <- paste0("chr", names(genome)) # annoying fix to get it to plot
seqTrack <- SequenceTrack(genome, genome = "hg38", chromosome = chrom,
                          noLetters = TRUE, fontcolor = c(`A` = "lightblue",
                                                          `G` = "#909AFD",
                                                          `C` = "#DAA400",
                                                          `T` = "#CF4918",
                                                          `N` = "gray"))

# bring in the bams
library(GenomicAlignments)

# set the constant scale
cov_scale <- c(0,40)

# make a dummy range for bulk, storm, and ss3xpress as there is no
# coverage here, but can't get the ylimits to plot when coverage is 0
dummy_range <- GRanges(seqnames = "chr6",
                       ranges = IRanges(start = from, end = to),
                       strand = "+")

# bulk
bulk <- readGAlignmentPairs("merged_bulk_100k_mapqfilt_nonannot_reads.bam",
                            strandMode = 2)
seqlevelsStyle(bulk) <- "UCSC"
bulk <- keepStandardChromosomes(bulk, pruning.mode = "coarse")
bulk_gr <- granges(bulk)
bulk_gr <- c(dummy_range,
             bulk_gr)
bulk_track <- AlignmentsTrack(bulk_gr, genome = "hg38",
                              chrom = "chr6",
                              name = "Bulk Total RNA-seq",
                              fill = "goldenrod", type = "coverage",
                              showTitle = FALSE,
                              stackHeight = 0.5)
displayPars(bulk_track) <- list(ylim = cov_scale)

# storm
storm <- readGAlignmentPairs("merged_storm_100k_nonannot_reads.bam",
                             strandMode = 2)
seqlevelsStyle(storm) <- "UCSC"
storm <- keepStandardChromosomes(storm, pruning.mode = "coarse")
storm_gr <- granges(storm)
storm_gr <- c(dummy_range,
              storm_gr)
storm_track <- AlignmentsTrack(storm_gr, genome = "hg38",
                              chrom = "chr6",
                              name = "STORM-seq", type = "coverage",
                              fill = "#63197FFF", showTitle = FALSE,
                              stackHeight = 0.5)
displayPars(storm_track) <- list(ylim = cov_scale)

# ss3
ss3 <- readGAlignmentPairs("merged_ss3_nonannot_reads_mapqfilt_100k.bam",
                             strandMode = 1)
seqlevelsStyle(ss3) <- "UCSC"
ss3 <- keepStandardChromosomes(ss3, pruning.mode = "coarse")
ss3_gr <- granges(ss3)
ss3_gr <- c(dummy_range,
            ss3_gr)
ss3_track <- AlignmentsTrack(storm_gr, genome = "hg38",
                             chrom = "chr6",
                             name = "Smart-seq3xpress", type = "coverage",
                             fill = "#8CCC98", showTitle = FALSE,
                             stackHeight = 0.5)
displayPars(ss3_track) <- list(ylim = cov_scale)

# vasa
vasa <- readGAlignments("merged_vasa_nonannot_100k_mapqfilt_reads.bam")
seqlevelsStyle(vasa) <- "UCSC"
vasa <- keepStandardChromosomes(vasa, pruning.mode = "coarse")
strand(vasa) <- "*"
vasa <- c(dummy_range,
          vasa)
vasa_track <- AlignmentsTrack(granges(vasa), genome = "hg38",
                              chrom = "chr6",
                              name = "VASA-seq", type = "coverage",
                              fill = "#CF5917FF", showTitle = FALSE,
                              stackHeight = 0.5)
displayPars(vasa_track) <- list(ylim = cov_scale)

# sstotal umi
sstotal_umi <- readGAlignments("merged_sstotal_umi_nonannot_reads_100k_mapqfilt.sorted.bam")
seqlevelsStyle(sstotal_umi) <- "UCSC"
sstotal_umi <- keepStandardChromosomes(sstotal_umi, pruning.mode = "coarse")
strand(sstotal_umi) <- "*"
sstotal_umi <- c(dummy_range,
                 sstotal_umi)
sstotal_track <- AlignmentsTrack(granges(sstotal_umi), genome = "hg38",
                              chrom = "chr6",
                              name = "Smart-seq-total", type = "coverage",
                              fill = "#8FC2FD", showTitle = FALSE,
                              stackHeight = 0.5)
displayPars(sstotal_track) <- list(ylim = cov_scale)

# bring in the polyA track bed file
polyA_track <- rtracklayer::import("hg38_polyA_min6_c2.bed.gz")
seqlevelsStyle(polyA_track) <- "UCSC"
polyA_track <- keepStandardChromosomes(polyA_track,
                                       pruning.mode = "coarse")
polyA_track$run_length <- unlist(lapply(polyA_track$name, function(x) {
  splt <- strsplit(x, "_")
  return(splt[[1]][3])
}))

polyA_track.roi <- subsetByOverlaps(polyA_track,
                                    roi)
polyaTrack <- AnnotationTrack(polyA_track.roi, name = "polyA-track",
                          shape = "box", from = from, to = to,
                          fill = "lightblue", showTitle = FALSE)

# build up the tracks to plot
plotTracks(list(ideoTrack, gat,
                bulk_track, storm_track, ss3_track, vasa_track, sstotal_track,
                polyaTrack,
                seqTrack, biomTrack), from = from, to = to,
           col.axis = "black",
           col.title = "black", background.title = "white",
           cex.axis = 1, font.size = 18,
           sizes = c(0.3, 0.3,
                     0.5, 0.5, 0.5, 0.5, 0.5,
                     0.2,
                     0.1, 0.3))
