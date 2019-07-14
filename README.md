# LEXY Toolkit

## Summary
This small collection of processing scripts was created for image analysis and quantitation of repurposed versions of the light-inducible nuclear export system (LEXY) optogenetic toolkit (Niopek D, Wehler P, Roensch J, Eils R, Ventura B Di (2016) **Optogenetic control of nuclear protein export.** _Nat Commun_, [doi:10.1038/ncomms10624](https://www.nature.com/articles/ncomms10624))

A summary of how to use these scripts in order is provided in the [workflow](workflow.md) file, including Jython functionality for inital processing of Leica LIF files in Fiji (ImageJ) and python scripts for quantitation. The provided test data includes a sample raw data file, and initial inputs for the python processing steps. Working copies of these scripts have also been provided as jupyter notebooks, which can be accessed via the [Binder](https://gke.mybinder.org/) badge above.

### ImageJ/Jython functionality
Briefly, the bounding region of interest (ROI) for all cells in a field of view is assigned on brightfield images using the [GDSC](http://www.sussex.ac.uk/gdsc/intranet/microscopy/UserSupport/AnalysisProtocol/imagej/toolsets) Cell Outliner plugin. Nuclei ROIs are assigned via automatic thresholding of the Hoechst image at each timepoint. Pixel-wise intensity information for both the cells and nuclei was then exported from all fluorescent channels for further analysis.

### Python functionality
Transfected cells are selected via thresholding the mean fluorescence intensity per cell and the corresponding nuclei automatically assigned via a maximum Euclidean distance threshold. These nuclei were then tracked over time via similar distance thresholdis. The mean fluorescence intensity of each nuclei 'track' over time is then exported to a single excel file for later plotting and analysis in graphical plotting software.

Additional functionality includes the filename_cleanup utility, before_after functions to compare the mean nuclear fluorescence before and after activation and the single_cell_plotter to visualise the mean nuclei fluorescence intensity over time for individual cells. Cells with aggregates were defined as those containing > 15% saturated pixels at time 0, and marked to allow removal from subsequent analysis.

This work will be prepared for publication, and additional details included here when available.

## Prerequisites

### Fiji/ImageJ
Use of this toolkit requires [Fiji](https://fiji.sc/) equiped with the [Bio-Formats](https://imagej.net/Bio-Formats) plugin and [GDSC](http://www.sussex.ac.uk/gdsc/intranet/microscopy/UserSupport/AnalysisProtocol/imagej/toolsets) Cell Outliner plugin. For additional details on installing ImageJ plugins, refer to the relevant [documentation](https://imagej.net/Installing_3rd_party_plugins).

### Python
The python processing scripts assumes a standard installation of Python 3.7. For specific package requirements, see the [requirements.txt](requirements.txt) file.


## Authors

* **Dezerae Cox** - *Initial work* - [dezeraecox](https://github.com/dezeraecox)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments and additional citations

* Templates for this package were adapted from [PurpleBooth](https://github.com/PurpleBooth)
* Kudos to the [2018 Asia-Pacific Advanced Scientific Programming in Python (#ASPP) Summer School](https://www.melbournebioinformatics.org.au/aspp-asia-pacific/) for giving me the condfidence and tools to tackle my first python package!
* ImageJ: Schneider, C. A.; Rasband, W. S. & Eliceiri, K. W. (2012), "NIH Image to ImageJ: 25 years of image analysis", Nature methods 9(7): 671-675, PMID 22930834 (on Google Scholar).
* Fiji: Schindelin, J.; Arganda-Carreras, I. & Frise, E. et al. (2012), "Fiji: an open-source platform for biological-image analysis", Nature methods 9(7): 676-682, PMID 22743772, doi:10.1038/nmeth.2019 (on Google Scholar).
