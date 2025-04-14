# STORM-seq inter-gene UMI collision analysis

Previous work has demonstrated that the occurrence of intragene UMI collisions is
near zero (https://www.nature.com/articles/s41587-021-00870-2), while inter-gene collisions
can be directly computed from the data. In this vein, we explored the rates of inter-gene
collisions across technologies. Using a simulation approach, we estimated the empirical mean
(mu ~ 6) and dispersion (theta ~ 6) parameters across genes as detailed in the methods of the paper. 
From there, we simulated inter-gene UMI collisions across fixed read depths using UMI sampling
with replacement to mimic sampling from the initial UMI pool during library prep. Given the approximate
gene detection rate of around 8k genes/cell across STORM-seq, VASA-seq, and SS3x at 100k reads/cell
we elected to use this gene detection rate as the approximate expected value for the obs/exp
calculation. Simulation results suggest the expected saturation as UMI length increases from
6 x N's in VASA-seq to the 16 x N's in Smart-seq-total. Inter-gene UMI collisions were defined
using primary alignments (not consider secondary or supplemental alignments) that map to
known/annotated genes in the Ensemble 101 gene annotations using STARsolo.

<img title="UMI collision detection" alt="UMI collision detection" src="/images/umi_collision_detection.png">
