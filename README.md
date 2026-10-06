# python_Voltage_Islet_Analysis
This repository includes the current Voltage Islet Analysis pipeline. This code is still in very active development and changes nearly every day. Only some of those changes will be commited to the GitHub if they prove to be substantial enough to justify a commit. It is also already used regularly for analysis of my own data. To maintain reproducibility the version that the data has been analyzed with is being documented.
# About AI usage
AI has been used for synthax and Debugging help. This is due to the fact that I am learning python while writing this analysis. Therefore I use a mixture of StackOverFlow and AI to find out how to phrase the needed commands.

# Intented workflow
The scripts are structured so that you are supposed to start with a tiff file. For now I resave my .czi and .orf into .tiff to feed them into the Pipeline.
For most applications (single Voltage/Calcium imaging) one does simply call the first script "Napari_ROI_picker.py". I prefer doing this in a conda environment because the last time when I didn't I destroyed my entire local python...which was a very sad day indeed. I repaired it and now I have learnt environments
