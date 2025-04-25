# STORM-seq manuscript analysis code and resources

The code and resources used in the STORM-seq manuscript have been deposited here.

## Important note about the cell/well names used

Within some of the provided scripts/manifest files/resource files the well names
may deviate from the downloaded GEO/SRA/ENA FASTQ files. The reason for this was
a miscommunication about the order of the UDIs. The files uploaded to GEO/SRA/ENA
have been adjusted to reflect the "correct" wells, will not need to go through any
of the "remapping" done in some of these scripts. For posterity, the remapping file
from the old wells (column 1) to the new/correct wells (column 2) of the `well_map.txt`.
