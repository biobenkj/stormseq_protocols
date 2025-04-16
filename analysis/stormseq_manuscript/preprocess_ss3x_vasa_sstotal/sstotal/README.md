# Smart-seq-total with UMIs pre-processing

Smart-seq-total with UMIs HEK293T data were downloaded from GSE151334 and pre-processed as previously described in https://github.com/aisakova/smart-seq-total with minor modifications. Briefly, read 2 reads were trimmed using default parameters with the following modifications: 1) `-j 4`, 2) `-u 6`, 3) `-a AAAAAAAAAA`, and 4) `-m 18`. Next, UMI containing reads (read 1) ending in a poly(T) sequence (‘TTT’) were extracted, and matched read 2 FASTQ files were re-paired using seqkit pair (v2.7.0). This was based on the recommendation from the Smart-seq-total authors:

```
For precise molecule count, one should filter R1 files (and R2 respectively) and keep only UMI-containing pairs (based on the presence of TTTTTTT at the end of each R1) (type 1 read above).
```

The extraction of reads in read 1 that end with the polyT motif was done using `filter_polyT.py`. 