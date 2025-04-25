# Demultiplex STORM-seq data from the manuscript

The data have been re-multiplexed for upload to GEO/SRA/ENA and live at: `GSE181544`.
Normally we wouldn't use the data as a re-multiplexed format as we have uploaded it, 
and instead use it on a per-cell basis, with the corresponding pair of fastq files - 
so something like `A1_L001_R1_001.fastq.gz A1_L001_R2_001.fastq.gz`. But in order to
minimize the number of files for upload, we have merged/re-multiplexed them by cell 
type and added a synthetic cell barcode with synthbar (https://github.com/jamorrison/synthbar).
The reality is that we took the i5 and i7 and concatenated them, similar to Smart-seq3xpress.

## Cell-type demultiplexing

Data can be downloaded from the above GEO/SRA/ENA location and converted back to fastq 
format (if you downloaded `.sra` format using sra-tools: https://github.com/ncbi/sra-tools).
As such, K-562, HEK293T, and RMG-2 fastqs contain multiple cells/wells within these merged 
fastq files. You can use the provided demultiplexing script (`demux_merged_storm.py`) in the 
following manner:

```
# example for K-562
python demux_merged_storm.py \
    --annot storm_annotated_metadata_with_md5.tsv \
    --read1 stormseq_K562_remux_R1.fastq.gz \
    --read2 stormseq_K562_remux_R2.fastq.gz \
    --cell-type K562 \
    --unassigned
``` 

## Checking file integrity

All demultiplexed files can then be checked against the MD5 checksum using `md5sum` with
the corresponding `storm_annotated_metadata_with_md5.tsv`. If they deviate, then something
may have happened during file conversion. If you run into this, please open an issue and let
us know. These MD5 checksums are on the _uncompressed_ FASTQ files to ensure byte-for-byte comparison
to the original files in `storm_annotated_metadata_with_md5.tsv` with those from GEO/SRA/ENA.

## Important note

Please trim off the first 16 bp from read 1 and read 2 in order to recover
the original/expected output from a STORM-seq run. Alternatively, you can skip
the demultiplexing altogether and specify the read architecture to include the
first 16 bp of both read 1 and read 2 to contain the relevant cell barcode.
