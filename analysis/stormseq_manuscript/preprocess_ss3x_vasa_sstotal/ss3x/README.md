# Pre-processing of SS3xpress data

## Demux and get sequencing depth of SS3xpress data
Generated fastqs from the unmapped bams using `ss3xpress_bam_to_fastq.py`. This produced 4 fastq files,
where I1.fastq.gz and I2.fastq.gz are the cell barcodes from the BC tag and if the UB tag is populated,
append that sequence with "perfect" base qualities (as they are lost in the unmapped bam) and reconstruct
the relevant 5' UMI and internal reads.

Generated per cell fastqs using the `ss3xpress_bam_to_sc_barcode_fastq.py`.
This script generates 4 fastqs per barcode from a known whitelist
Please note that there are no correction features for barcode extraction,
this needs to be done prior to execution.

After per cell fastq generation, cells were calculated for read depth using awk.
As an example:
First calculate read depth using any one of the 4 fastqs per barcode - awk 'END { print NR / 4 }' filename.txt
Then go through and find files/cells with at least 100k reads - awk '$2 > 100000 { print $1 }' filename.txt
Column 1 here is the file name/cell barcode and column 2 is the read depth calculated above.
Note: dropped 2 cells from the expected 96.

Subsampling was done with seqtk for each of the 4 fastq files setting a seed to 100 for
seqtk sample. These subsampled files were written directly to the appropriate directory.

```
# something like
seqtk sample -s100 fq1.fq 100000 | gzip -c > fq1.100k.fq.gz  
```

## Separate 5' UMI and Internal read fragments

To separate the 5' UMI and internal read fragments, the `sep_ss3xpress_read_types.py` script was used
to search for the conserved `ATTGCGCAATG` sequence to delineate the 5' UMI reads, similar to zUMIs
`fqfilter_v2.pl` script. Reads not possessing this conserved sequence were considered to be internal
reads. 