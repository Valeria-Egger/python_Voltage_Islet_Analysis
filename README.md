# python_Voltage_Islet_Analysis
This repository includes the current Voltage Islet Analysis pipeline. This code is still in very active development and changes nearly every day. Only some of those changes will be commited to the GitHub if they prove to be substantial enough to justify a commit. It is also already used regularly for analysis of my own data. To maintain reproducibility the version that the data has been analyzed with is being documented.

This code has been developed by me to analyse experimental data. Therefore the reasoning of every step will be explained in this ReadMe. It also means that the reasoning of every step is still very much under constant interrogation and may change any time. Since I started learning python 6 months ago, some of the concepts and design choices may reflect on this learning curve. \
**About AI usage**\
AI has been used for synthax and Debugging help. This is due to the fact that I am learning python while writing this analysis. Therefore I use a mixture of StackOverFlow and AI to find out how to phrase the needed commands.

## Intented workflow
All concepts described here will be discussed in more detail below.
It is recommended to run the analysis in a dedicated environment installed with python 3.10. Every other python version has not yet been systematically tested and might lead to dependency related problems. The dependencies will be documented once the scripts have been consolidated accordingly.
### Napari ROI picker
One starts by running the file "Napari_ROI_picker.py". This expects a .tiff timeseries with the dimensions HxWxT. It is not currently possible to load a 4D image (HxWxZxT) in this script only one plane for selection is supported. The script provides an interactive plot to quickly check traces already at this step for potential activity but also to check Z-dirfts or motion artefacts and potentially correct them. For motion correction a Matlab script is used for now. 
The script saves all ROIs in a .npy format which is then used for further Analysis.
### Load data
Once the .npy files are obtained one proceeds with the file "load_data.py". This script includes normalization of the trace, dff calculations as well as Z-scores that will later be used for peak detection. Furthermore this is the place to initialize variables for perfusion setups, frame rates and conditions. Also it is the place to tell the script whether one uses Calcium or Voltage data to analyze which changes the calculations of the aforementioned metrics. The data is stored in a dataframe with the fields "dataset" (for the name of the dataset, can be set by the user), "trace_id" (to distinguish individual traces), "condition" (for example beta-ablated and control, can be set by the user), "fs" (framerate, can be set by the user), "time" (depends on framerate), "signal" (the raw traces), "background" (the raw background traces), "baseline" (the calculated baseline of each trace), "dff" (the calculated dff of each trace), "z" (the calculated z-scores of each trace), "glucose" (the perfusion conditions for each timepoint of each trace)
This dataframe is stored as a pickle for further processing. 
### detect events and calculate summary metrics
After having the first normalization calculations one may proceed with the file "detect_events.py". This file expects the aforementioned pickle file format for the dataframe and tries to find events based on it. For event detection of the Calcium data the user has to manually set enter and exit threshold as well as the gap between events. For most of the voltage data I have, an enter threshold of Z-score -3 is sufficient to detect most downward deflections. The exit threshold of -2 makes sure that several deflections that may come from the same event stay merged as one event. The gap makes sure that even if the exit threshold should be too high, too small of a gap will not result in a new event. 
These metrics highly depend on the form of the data one has and may be adjusted manually. I am aware that this form of detection holds many risks for both false positives and negatives, however I will adress this question in the future. For now it works well enough when combined with human observer.
The peak detection for Calcium is different than for Voltage imaging therefore this separation exists in this script as well. It also depends on the framerate set by the user with the logic that everything over 100Hz is considered voltage everything under it is considered Calcium. This will also be changed in the future.
The event dataframe includes the following fields: "dataset", "condition", "trace_id", "fs" (taken from the other dataframe, fs is the framerate), "event_idx" (Index for each Event) "start_idx", "end_idx" (Indices where an event starts or ends), "duration seconds" (duration in seconds), "duration total" (duration without a formal time metric assigned), "amplitude" (amplitude of each event), "event area" (cumulative sum (AUC) of just the frames marked as an event), "glucose" (perfusion dependent glucose concentration for each event).
After this the summary dataframe is computed with some more metrics, the dataframe includes the following fields:
"dataset", "condition", "trace_id", "fs" (metrics taken from the other dataframes), "n_events" (count of events per trace), "mean amplitude", "median amplitude" (mean and median of amplitude per trace), "max amplitude" (maximal amplitude in that trace), "mean duration", "median duration" (mean and median of the duration seconds variable (so the time metric assigned duration)), "total area" (sum of all event areas), "event frequency" (frequency of oscillatory components of events using a Morlets wavelet transform), "dominant frequency" (dominant frequency of the entire dataset, using Morlet wavelet transform), "dominance ration", "frequency z-scores" (significance metrics for the reliability of the dominant frequency calculation), "glucose" (perfusion dependent glucose concentration at each timepoint), "AUC" (area under the curve for the entire trace not just events for each differing glucose concentration), "Event area" (sign independent area under the curve calculations for the entire trace for each differing glucose concentration).
This dataframe is being saved as a pickle file and can be reused for plotting in another file. The plots in this file mainly serve debugging and exploration purposes. 

## Information about Calculations
Aquisition parameters:
The images obtained were imaged with a frequency between 55 Hz and 780 Hz. However most images are either imaged with 100, 200 or 300 Hz respectively, which is the main focus of later processing steps. 
The reasoning for this range of frequencies is mainly taken from other papers and trial and error (hence the huge range of frequencies). However to formalize this step I calculated the Nyquist limit for my endocrine voltage events. Literature claims that endocrine voltage events are between 50 and 200 ms. \
```math
f_s=\frac{2}{T} =\frac{2}{0.05} = 40
```
```math
f_s=\frac{2}{T}=\frac{2}{0.2}=10
```
Since we like an oversampling of 3-5 times the Nyquist limit to actually resolve events that would translate to 120-200 Hz (for 50 ms events). In case that there are shorter events higher frequencies were also applied. Overall the band of 100-300 Hz should resolve the expected events.
Most of the time 4x4 Binning was applied to the images to boost signal to noise ratio already in the acquisition step. In some cases a 2x2 Binning was deemed as sufficient to further process the timeseries. This may interfere with cell identification.

**Quality control:**\
A small control step was implemented to check pictures for motion, focus drift, expression and bleaching.\
**Motion score**:\
For each frame the absolute difference between the pixels in this frame and the previous frame is calculated and the averaged over all pixels to collapse the absolute difference into one number per timepoint. The absolute difference has the advantage that positive and negative signs will be ignored and cannot cancel themselves out. This should mainly detect big changes that affect many pixels in the frame, rather than intensity changes that my signals would produce (but only in a subset of pixels that should be mostly cancelled out by averaging over all pixels including background that should not move at all). Due to noise and small fluctuations the motion score will never be zero. Therefore the mean of all timepoints is normalized to the mean expression level in the timeseries.
Source:
(https://opencv.org/autofocus-using-opencv-a-comparative-study-of-focus-measures-for-sharpness-assessment/#h-explanation-of-different-focus-measurement-techniques) \
**Focus score:**\
For focus score the Laplacian variance is being used as a metric on how blurry the image is. It measures the variance of the Laplacian response, assuming that sharper edges produce higher variance while blurry edges produce lower variance.

```math
\sum_{x,y} |\nabla^2 I (x, y)
```
with $\nabla^2(x, y)$ is the Laplacian \
Source:
(https://opencv.org/autofocus-using-opencv-a-comparative-study-of-focus-measures-for-sharpness-assessment/#h-explanation-of-different-focus-measurement-techniques) \
**Expression Score:**\
Simply calculates the pixel intensities of the whole image, if the image is brighter there is bright signal. Also the expression level of the image is being used as a normalization for the other metrics. \
**Bleaching rate:**\
Bleaching is estimated as a linear regression over the intensity of the pixels. The slope of this fitted linear regression is then stored as the bleaching rate. I use a linear regression here rather than a bi-exponential fit (which according to some sources should estimate bleaching better) because some of my traces fail on this fit, while the linear regression seems to be more stable for me at least. While I can measure bleaching and correct for it in the trace, it is impossible to recover worsening SNR ratio introduced by this bleaching, limiting my acquisition time to around 5-10 minutes (depending on signal strength, microscope that I use and other things). 

```math
I_t = \alpha + \beta t + \epsilon _t
```
with $\alpha$ = intercept, $\beta$ = slope (bleaching rate), $\epsilon _t$ = residuals \

After passing these metrics the dataset is further processed. \
**Motion correction:** \
Motion correction:
As motion correction NormCorre is used. NormCorre is a rigid and non-rigid motion correction algorithm that splits the image into overlapping patches and arranges each patch to a template it creates out of the first few hundred frames of the time series. It has been developed for in vivo motion correction and can therefore correct for breathing, heartbeat and blood flow but also deformations of the islet in the ex vivo culture. \
Source:
(https://www.sciencedirect.com/science/article/pii/S0165027017302753) \
**Preprocessing steps (Voltage imaging):** \
The mean intensity of each ROI is stored as an 3D array with form HxWxT.
For baseline estimation a lowpass filter is being applied to the trace. This filters out all frequencies above a certain cutoff that is dependent on the sampling frequency in the following manner: \
Cutoff: $\frac{0.05}{\frac{f_s}{2}}$ \
To be totally honest here the 0.05 was simply chosen by testing several numbers that were mentioned by the internet for “useful” lowpass filtering ranges and taking the best one. This needs adjustment for sure. 
After getting the baseline now only consisting of the frequencies below the threshold the trace is normalized using the baseline.\
```math
dff = (\frac{F_t - F_0}{F_0})
```
with dff = normalized trace, $F_t$ = trace at time t, $F_0$ = baseline \

This normalization step already helps with bleaching, removes artefacts and let’s the trace start at 0. The changes are easily converted in percentages by multiplying them with 100, which makes comparisons easier and amplitude more intuitive. 
According to literature voltage events from this sensor (Voltron2) are most commonly found in the 3-8% range, sometimes they can be higher (up to 15-20%) but the expected range is mostly 3-8%. \
**Peak Detection (Voltage Imaging):** \
For peak detection a Z-score algorithm was implemented. 
First of all the median of the trace is added to the dff to restore the baseline of the trace and prevents events from contaminating the Z-scores.
Then a quiet window is defined using the upper 50 percent for the residual trace that is the quiet window. This is because voltron2 produces downwards spikes therefore there is a higher probability that anything in these upper 50% of the trace are not signals but noise.
This could however be risky if we look at repolarisation, however I have yet to find a better definition for a quiet window.
Standard deviation was tried for the traces but those produced huge Z-scores for silent windows and more balanced but smaller Z-scores for windows with actual events. Therefore the Mean absolute deviation was used.
```math
MAD = \frac{1}{n} \sum_{i=1}^n |x_i - \mu |
```
with $x_i$ = each datapoint, $\mu$ = mean of datasets, $n$ = number of observations, $|x_i - \mu |$ = absolute deviation of each datapoint from the mean. \
this assumes a gaussian distribution of the datapoints that I have not yet in fact shown to be actually there to be honest. However apparently many people use this exact method. \
To then estimate the standard deviation from the MAD (with a gaussian distribution) the following formula was used: 
```math
\sigma = \frac{MAD}{0.6745} = 1.4826 \cdot MAD
```
This estimation for the standard deviation allowed for a more stable Z-score calculation later on. However I still have some problems with rather high Z-scores in silent windows that may inflate the positive rates in these windows. I have to yet find a solution for this problem. One possibility might be to couple it to a hard threshold that has to be met before calling a spike an event. \
To calculate the Z-scores the typical formula was used: \
```math
Z = \frac{(x_i - \mu)}{\sigma}
```
with $\mu$ = mean of the dataset, $x_i$ = datapoint at each timepoint, $\sigma$ = standard deviation estimate using the MAD \
Right now in order to call something an event hard threshold are defined. For onset of an event a threshold of -3 has to be met. For this event to stop being an event the signal has to return to a Z-score of under -2. Whenever the trace drops below the -3 mark the detector starts calling $x_i$ an event. This event stops only when the trace returns to over -2. It has been discovered that the different threshold for onset and offset help with the small fluctuations within the trace. If peaks are very close to the threshold one consecutive peak might be split in several small downward spikes without a defined offset score that is more lenient than the onset score. If peaks are being found in close proximity of 50 timepoints of each other they are being merged. This should prevent again splitting of peaks and resembles physiology where a cell cannot depolarize shortly after a depolarization. \
It should be noted that this peak detection is still very rudimentary, since it depends on hard threshold only, which makes it’s adaptability limited. I have yet to find another way to do this. \
From this peak detection events are being defined and fed into an event dataframe which gives each event in a trace an index and calculates characteristics of these events. \
**Preprocessing Calcium Imaging:** \
For Calcium imaging a simple min-max scaling is enough for normalization and dff calculations:
```math
dff = frac{(F_t-min)}{(max-min}
```
with min = Minimum value in the trace, max = maximum value in the trace, $F_t$ = trace at each timepoint \
The Z-trace is calculated the same way as in the Voltage version:
```math
Z = \frac{(x_i - \mu)}{\sigma}
```
with $\mu$ = mean of the dataset, $x_i$ = datapoint at each timepoint, $\sigma$ = standard deviation \
**Peak Detection Calcium Imaging:** \
For Calcium imaging the peak Detection is taken mainly from the following post from StackOverflow: (https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/43512887#43512887)
**Dataframe calculations:** \
Following characteristics are being defined on the event level: \
Duration/total Duration: $d = (end - start) $ \
Duration in seconds: $d_s = \frac{(end - start}{f_s}$ \
with $f_s$ = sampling rate \
Amplitude : Minimum of the trace (negative going indicator) for the event at that timepoint. \
Area: $area = \sum x_i$ \
with $x_i$ = datapoint at each timepoint \
cumulative sum for Area under the curve estimation. \
Furthermore the event dataframe saves start and end indexes of each event, trace ID, event ID and sampling rate.
Then to make comparisons not only between events but also between cells (traces) and whole datasets a summary dataframe has been implemented measuring the following characteristics of the traces: \
Count of events \
Mean/median amplitude \
Mean/median duration \
Total area under the curve \
Event frequency: $Event_freq = \frac{count of events}{total time}$ \
Dominant frequency: computed using the cwt command in python that computes a Morlet wavelet transform (https://www.sciencedirect.com/science/article/pii/S0888327007000994) \
This lets one extract the wavelet coefficients of the trace, which can then be converted to power using the following calculation: $power = |c^2|$ \
which represents the energy present. Then the frequencies are being capped to only show relevant frequencies rather than noise that might dominate the imaging. Then the mean power per scale is calculated and from that the maximum shows the dominant scale or frequency of the entire dataset. \
Dominance ratio for that frequency: $dominance ratio = \frac{P(f_{dom})}{mean(P(f))}$ \
Z_scores of that frequency: $Z = \frac{p - \mu f}{\sigma f}$ \
with $\mu f$ = mean of frequencies, $\sigma f$ = standard deviation of frequencies
These dataframes are being saved as pickles (.pkl) to preserve python objects and avoid loading the entire timeseries into the storage space when I can simply work on the dataframe itself.
Furthermore more characteristics may be calculated once I find the need for them. 
 ## Limitations and future plans
 -hardcoded saving variables will be removed/exchanged with an UI \
 -hardcoded entry paths may be considered to be exchanged with an UI \
 -hardcoded thresholds and variables will be reconsidered \
 -peak detection on Calcium needs refinement \
   -there is a bug in the peak detection for Calcium overall making the output of n_events always 0 although peaks are being detected \
  -peak detection on Voltage traces is still crude and may be improved \
  -dual imaging napari ROI picker will be improved \
  -plotting and analysis will be separated \
  -more tests may be carried out to see performance also on Calcium datasets \
  -Dependencies will be evaluated \
  -...
