# invitRo

Prioritization of LC-HRAM-MS non-target screening features acquired with a pre-defined human liver microsomes 
in vitro metabolism assay as described by Mortelé O., Vervliet P. et al (https://doi.org/10.1016/j.jpba.2018.02.032).
An overview of the samples in the assay can be found in the supplied CSV file (HLM_samples.csv)
A visual representation can be found in the supplied image (HLMassay.png)

The purpose of the script is the prioritization of features previously picked from acquired LC-HRAM-MS. 
In our workflow we converted instrument data files to mzXML data format using MSConvert (ProteoWizard, 
http://proteowizard.sourceforge.net/tools.shtml) and imported them to MZMine (http://mzmine.github.io/).

First, m/z features were detected using the centroid algorithm followed by a chromatogram building step. 
Resulting chromatograms were deconvoluted using the noise amplitude algorithm. 
Further simplification occurred by deisotoping, keeping the lowest m/z value as the representative isotope. 
During the deisotoping, chromatographic peaks were also filtered according to peak width: only peaks with a width 
between 0.05 and 1 min were retained. Next, peaks were aligned across samples using the random sample consensus (RANSAC) 
alignment algorithm. Finally, using the Same RT and m/z range gap filler algorithm any missing peaks were re-iteratively 
extracted. 

The final peak list was exported to a CSV file with a structure equal to the sample file availble in the repository. 
-> Sample-InVitroMassList R.csv

The script allows the user to create a PCA of the different samples and to prioritize features using a visual representation of each sample group in a volcano plot. A subset of the features of interest can be exported to a csv file for further identification. 

DOI: 10.5281/zenodo.1252144
