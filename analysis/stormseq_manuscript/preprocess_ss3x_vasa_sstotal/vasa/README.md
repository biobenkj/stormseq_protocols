# VASA-seq pre-processing

Passage matched K-562 and HEK293T VASA-seq data were generated to directly compare to STORM-seq.
Libraries were prepared by Single Cell Discoveries (https://www.scdiscoveries.com/services/vasa-seq/) in an effort to avoid library preparation bias as we have not generated VASA-seq libraries internally. Reads were sequenced on a NextSeq using a 1x75 bp kit and two plates in total were generated. Quality control metrics received from Single Cell Discovery are included for convenience. Plates are evenly split between K-562 (first 192 wells per plate) and HEK293T (remaining wells to 384). Note that initial sequencing depth was not met, so libraries were resequenced and merged together (available on GEO). Negative controls are found in `neg_controls.txt` and should be filtered out for downstream analyses. Bioanalyzer traces and SCD generated QC metrics can be found in the `vasa_qc` directory. In total across both 384-well plates, >92% of single cells had libraries passing QC thresholds with >1000 endogeneous reads, and an average raw read depth of 179,500 reads/cell.

## Demultiplexing

VASA-seq data were demultiplexed to single cell FASTQ files using provided scripts from https://github.com/hemberg-lab/VASAseq_2022 and specific scripts can be found in `1_demux`. 

## Trimming

VASA-seq data were pre-processed and trimmed using cutadapt (v4.4) through provided scripts at https://github.com/hemberg-lab/VASAseq_2022 and can be found in `2_trim`. After trimming, the UMI and cell barcode were appended to the 5’ end of the read using `fix_vasa_fqs.py`. Importantly, in the provided `concatenator.py`, a conversion error occurs when extracting and appending the UMI and cell barcode where a base quality score of 14 (ASCII character ‘/’) is unintentionally converted to a base quality score of 46 (ASCII character ‘O’). As a result, during the addition of the UMI and cell barcode to the 5’ end of the read, these base qualities were reverted to their original value.


<img title="VASA-seq bioA" alt="VASA-seq bioA" src="/images/vasa_bioa_cdna.png">
