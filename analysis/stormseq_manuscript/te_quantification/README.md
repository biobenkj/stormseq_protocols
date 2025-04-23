# Transposable Element (TE) quantification

The following workflow was modeled after Shao and Wang, 2020 Genome Research
(https://genome.cshlp.org/content/31/1/88). We deviate in the sense that we do
not attempt to perform either de novo or guided TE-transcript assembly, and restrict
to annotated TEs found within intergenic and intronic space (repeat masker and Ensembl 101).
The full pipeline used can be found (https://github.com/huishenlab/TE_quantification_pipeline).
If you prefer snakemake, we have adapted this pipeline to be used with snakemake (https://github.com/huishenlab/TE_quantification_snakemake).
After resource generation (or downloading the resources we used from zenodo), locus-level
TE expression was quantified using counts per million (CPM) and featureCounts, and
converted to an observed/expected value (see methods of the STORM-seq paper or the above github links).

## Filtering and plotting

Observed/expected values were plotted using ggplot and found in the 

<img title="TE exp HEK" alt="TE exp HEK" src="/images/te_quant_hek.png">
