# python_Voltage_Islet_Analysis
This repository includes the current Voltage Islet Analysis pipeline. This code is still in very active development and changes nearly every day. Only some of those changes will be commited to the GitHub if they prove to be substantial enough to justify a commit. It is also already used regularly for analysis of my own data. To maintain reproducibility the version that the data has been analyzed with is being documented.

This code has been developed by me to analyse experimental data. Therefore the reasoning of every step will be explained in this ReadMe. It also means that the reasoning of every step is still very much under constant interrogation and may change any time. Since I started learning python 6 months ago, some of the concepts and design choices may reflect on this learning curve. 
## About AI usage
AI has been used for synthax and Debugging help. This is due to the fact that I am learning python while writing this analysis. Therefore I use a mixture of StackOverFlow and AI to find out how to phrase the needed commands.

## Intented workflow
### Napari ROI selection
All concepts described here will be discussed in more detail below.
It is recommended to run the analysis in a dedicated environment installed with python 3.10. Every other python version has not yet been systematically tested and might lead to dependency related problems. The dependencies will be documented once the scripts have been consolidated accordingly.
One starts with running the script "Napari_ROI_picker.py" which is sufficient for all single Voltage/Calcium analysis to extract traces out of. If you already have extracted traces in a .npy format you are welcomed to skip this step. 
For dual imaging of Calcium and Voltage using the script "dual_napari_picker.py" helps in providing interactive plots for both signals and separates the signals by itself. This might help in quickly checking of whether or not your dual imaging produces expected colocalization of events. However it is also perfectly valid to run both traces through the normal Napari-picker and load them separately since the loading itself does have to occur separately for these traces. 
Please also note that the background for dual imaging is for both traces the same. 
