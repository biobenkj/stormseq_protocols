# Using Existing QuPath Classifiers

An existing QuPath project with images and detections (not cell objects) is required to use the object classifiers. The objects additionally need to have the same measurements with the same channel names. Alternatively, edit channel names in classifier json to match channel names in project.

## Required measurements

- CD31: Mean
- CD31: Median
- CD31: Min
- CD31: Max
- CD31: Std. Dev.
- CD45: Mean
- CD45: Median
- CD45: Min
- CD45: Max
- CD45: Std. Dev.
- PanCK: Mean
- PanCK: Median
- PanCK: Min
- PanCK: Max
- PanCK: Std. Dev.

## How to use

Copy the json files of the object classifiers (or whole folder) into the classifier folder of the QuPath project. The classifier can now be used through the GUI in `Classify > Object Classifier > Load Object Classifier` or through scripting using the command `runObjectClassifier("name of classifier")`.

Full script to create detections and run classifier can be found [here](/analysis/stormseq_manuscript/cycif/CellPose.groovy).

The CellPose python environment and QuPath extension are required for using the full script. Links to resources below:

[QuPath 6.0](https://qupath.github.io/)
[CellPose](https://github.com/MouseLand/cellpose)
[CellPose QuPath Extension](https://github.com/BIOP/qupath-extension-cellpose)

***Disclaimer***
*The object classifiers were specifically trained on the CellDIVE dataset included in the Johnson, B.K. et.al. publication. This model may not directly transfer to data collected from different microscopes, different staining protocols, or different imaging parameters for other samples stained for CD31, CD45, and/or PanCK. We recommend training your own object classifier for your own imaging data using the same approach described in the methods section.*
