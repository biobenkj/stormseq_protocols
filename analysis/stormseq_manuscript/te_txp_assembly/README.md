# TE-derived transcript stitching and quantification

Prior to any analyses, TEProf3 was used to generate the reference data files as
described in the TEProf3 documentation using Ensembl 101 annotations and the
latest repeat masker file (see manuscript methods). The `repeatmasker_sorted.txt`
file was modified to only include TE candidates that overlapped the nascent RNA 
(PRO-cap/GRO-cap) peaks (distal and proximal) from PINTS (https://pints.yulab.org/).
This modified file can be found in the `resources` directory.

## TEProf3 TE-derived transcript stitching calls

### Pseudobulk

The following was used to generate the `teprof3_output_TE_transcript_consensus.gtf`
file that was then used to quantify single-cell TE-derived transcripts.

```
teprof3 --manifest pb_sample_manifest.txt -am 1 -at 40 -fm 1 -ft 8 -tt 40 -ql 150 -ki
# -am 1 how to run transcript de novo assembly: 1 = short-read only
# -at 40 number of threads used for assembly: 40
# -fm 1 how to filter TE-derived transcripts: 1=only short read
# -ft 8 threads for filtering transcripts: 8
# -tt 40 number of threads for TACO: 40
# -ql 150 read length of libraries: 150
# -ki keep intermediate files
```

### Single-cell TE-derived transcript quantification

The following was used to generate the raw TE-derived transcript quantification
table (`teprof3_output_quantification.TE.tsv.gz` and associated annotations
`teprof3_output_filter_transcript_TE_transcript_consensus.tsv.gz`). An important
note is that originally the stringtie quantification call was hard coded for
unstranded quantification, while STORM-seq is reverse stranded. A TODO is to further
modify the codebase to allow for that to be an option instead of hard coding it. The
quantifications presented in the manuscript was modified to accommodate the reverse-stranded
read structure of STORM-seq.

```
teprof3 --manifest sc_sample_manifest.txt -ki --guided teprof3_output_TE_transcript_consensus.gtf -ql 150 -s 40
# --guided run TEProf3 in guided mode to use a consensus TE-derived transcript GTF for quantification
# -s 40 number of samples to be processed together at the same time: 40
```

### Quantification of TE-derived transcripts in bulk total RNA-seq

A similar call to the above STORM-seq TEProf3 guided mode was used.

```
teprof3 --manifest bulk_sample_manifest.txt -ki --guided teprof3_output_TE_transcript_consensus.gtf -ql 50
```

## Processing and plotting of K-562 TE-derived transcripts

Please see the `teprof3_postprocessing.R` script for details. Briefly, bulk and STORM-seq
data were filtered as recommended by TEProf3 and expressed in at least 10% of cells, respectively.

## LTR1A2-PURPL analysis and plotting

Single-cell BAMs out of STARsolo were overlapped with a merged GTF of canonical PURPL transcript isoforms
and the TEProf3 stitched LTR1A2-PURPL transcripts. BAMs were subset to reads overlapping this region with
samtools to increase computational speed. FeatureCounts was then used to tease apart reads that likely belong
to either the canonical isoforms or the TE-derived transcript (LTR1A2-PURPL). Following featureCounts, the
`parse_fc_te_and_known_txps.py` was used to redistribute reads into the canonical PURPL or TE isoforms.
Downstream processing was done with a slightly modified Gviz in R (https://github.com/biobenkj/Gviz) where
sashimi plots are generated as a proportion of the joint spanning reads. Plotting and analysis can be found
in the `purpl_plotting_gviz.R`.

