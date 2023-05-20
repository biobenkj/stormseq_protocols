# STORM-seq protocols for library preparation and analysis

This repository contains the current and archived STORM-seq protocols
for library construction and analysis. 

### Library construction
You can expect this process to take approximately
8-10 hours to go from single cells to prepared library. Please feel
free to open an issue if you encounter any sort of problem, have a
question, or have a tip/trick that you'd like to share/see implemented.

### Analysis
STORM-seq can ultimately be analyzed as if it were a "bulk" RNA-seq
protocol and is functional with tooling that supports this kind of 
input. Additionally, we recommend and support analysis using Kallisto|Bustools
and STARsolo. The required software and scripts to run these tools are
found in the 'analysis' folder/directory. Additionally, we are working
on producing docker and singularity containers for distribution as well.
As with the library prep documentation, please feel free to open an issue
if you encounter a problem or have suggestions for feature development.

*Thanks for trying out STORM-seq!*

#### Key kit components to purchase:

    * STORM-seq kit: Takara Bio cat: 634751
    * Takara Bio UDI sets A-D: Takara Bio cats: 634752, 634753, 634754, 634755
    * ERCC spike-in mix: Thermo Fisher cat: 4456740

#### Plate types to purchase:

    * Eppendorf twin-tec lobind skirted 384-well plates: Eppendorf cat: 0030129547
    * Eppendorf twin-tec lobind skirted 96-well plates: Eppendorf cat: 0030129512
    * If using SPT Mosquito as automation: SPT 384 LVST plate cat: 4150-05829

*Please refer to the latest protocol for remaining reagents needed*

#### STORM-seq workflow diagram

<img title="STORM-seq workflow" alt="STORM-seq workflow" src="/images/storm_workflow.png">
