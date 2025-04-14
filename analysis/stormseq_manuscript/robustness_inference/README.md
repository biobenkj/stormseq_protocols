# Transcript detection and variance analysis

In an effort to quantify transcript detection and variance across technologies, we developed
a simulation/resampling framework to attempt to begin to answer this question. First, we 
subsampled single cells using seqtk to the maximum read depth possible across technologies, and were
limited to 150k reads/cell in VASA-seq. After subsampling, to simulate 10 separate experiments using
these single-cell expression profiles, cells were subsampled with 10 random seeds (`rand_seeds.txt`)
to 50k reads/cell. These were then pseudo aligned and quantified using kallisto|bustools at the transcript
level as described in the methods of the paper. Per-cell bootstrapped quantification was then processed
and summarized as described in the methods of the paper and associated R script (`bootstrap_tech_txp_robustness_inference.R`).
Transcript analyses were carried out across commonly detected transcripts in all technologies. Transcripts
below 500 bp based on Ensemble 101 annotations were not used as these are typically removed during library prep.

<img title="Txp detection robustness" alt="Txp detection robustness" src="/images/txp_robustness.png">
