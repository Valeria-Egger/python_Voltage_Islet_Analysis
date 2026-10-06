import napari
import skimage.draw as draw
import numpy as np
import tifffile as tiff
import cv2
from skimage.morphology import closing, disk, remove_small_objects
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import savgol_filter
import pandas as pd
import math
import pywt
from scipy.optimize import curve_fit
from matplotlib.widgets import SpanSelector
from magicgui import magicgui
import seaborn as sns

#again no exploratory science in here

df_traces = pd.read_pickle("df_traces_plateau2.pkl")
enter_threshold = -3
exit_threshold = -2.5
min_length = 50
gap = 100

# Source - https://stackoverflow.com/a/43512887
# Posted by Roman Kiselev, modified by community. See post 'Timeline' for change history
# Retrieved 2026-09-11, License - CC BY-SA 4.0

#!/usr/bin/env python
# Implementation of algorithm from https://stackoverflow.com/a/22640362/6029703
import numpy as np
import pylab

def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))
#of course the peak detection for Calcium and Voltage is different
#so here is the one for Calcium
#adapted from the internet and my previous Matlab script
def Calcium_events(z_trace, lag=100, threshold = 3.0, influence = 0.1):
    out = thresholding_algo(z_trace, lag, threshold, influence)
    sig = out["signals"]
    sig = (sig==1).astype(int)
    starts = np.where(np.diff(sig)==1)[0]
    ends = np.where(np.diff(sig)==-1)[0]
    if len(starts)==0 or len(ends)==0:
        return []
    if sig[0]==1:
        starts = np.insert(starts,0,0)
    if sig[-1]==1:
        ends = np.append(ends, len(sig)-1)
    return list(zip(starts, ends))
#and here is the Voltage peak detection
#it is a lot less sophisticated in concept, just a MAD based Z-score thresholding
def detect_events(z_trace, fs, enter_threshold, exit_threshold, min_length, min_gap):
    in_event = False
    signals = np.zeros_like(z_trace)

    for i in range(len(z_trace)):
        if not in_event and z_trace[i] < enter_threshold:
            in_event = True
            signals[i] = 1
        elif in_event and z_trace[i] > exit_threshold:
            signals[i] = 0
            in_event = False
        else:
            signals[i] = signals[i-1]

    event_starts = np.where(np.diff(signals)==1)[0]
    event_ends = np.where(np.diff(signals)==-1)[0]
    if len(event_starts) == 0 or len(event_ends) == 0:
        return []
    if signals[0]==1:
        event_starts = np.insert(event_starts, 0, 0)
    if signals[-1] == 1:
        event_ends = np.append(event_ends, len(signals)-1)

    merged_events = []
    current_starts = event_starts[0]
    current_end = event_ends[0]

    for s, e in zip(event_starts[1:], event_ends[1:]):
        if s - current_end < min_gap:
            current_end = e
            continue
    
        merged_events.append((current_starts, current_end))
        current_starts = s
        current_end = e

    merged_events.append((current_starts, current_end))
    return merged_events

#at some point I really liked wavelets so of course they get to make an appearance
def coeffs_func(trace, scales, sampling_dt):
    coeffs, freqs = pywt.cwt(trace, scales, "morl", sampling_period = sampling_dt) 
    return coeffs, freqs
#and at some other point I realized you can even do stuff with wavelets so that got to make an appearance as well
def dominant_frequency(trace, scales, sampling_dt, scale_min, scale_max, fc = 0.8125):
    coeffs, freqs = coeffs_func(trace, scales, sampling_dt)
    power = np.abs(coeffs)**2
    scale_mask = (scales >= scale_min)&(scales <= scale_max)

    band_power = power[scale_mask].sum(axis = 0)
    mean_power_per_scale = power[scale_mask].mean(axis=1)

    index_dom = np.argmax(mean_power_per_scale)
    scale_dom = scales[index_dom]

    freqs_band = freqs[scale_mask]
    f_dom = freqs[index_dom]

    peak = mean_power_per_scale[index_dom]
    background_f = np.mean(np.delete(mean_power_per_scale, index_dom))
    dominance_ratio = peak/background_f

    mu_f = np.mean(mean_power_per_scale)
    sigma_f = np.std(mean_power_per_scale)
    z_score_f = (peak-mu_f)/sigma_f
    #here I actually compute 2 significance metrics because I wasn't exactly sure about just one 
    return f_dom, dominance_ratio, z_score_f


event_rows = []
for tid in df_traces.trace_id.unique():
    df_t = df_traces[df_traces.trace_id == tid]
    z = df_t["z"].values
    fs = df_t["fs"].iloc[0]
    dataset = df_t["dataset"].iloc[0]
    condition = df_t["condition"].iloc[0]
    #again here is the Voltage/Calcium split
    if df_t["fs"].iloc[0] > 99:
        events = detect_events(z, fs, enter_threshold, exit_threshold, min_length, gap)
    else:
        events = Calcium_events(z)

    for event_id, (start, end) in enumerate(events):
        total_duration = end-start
        duration_s = (end-start)/fs
        amp = df_t["dff"].iloc[start:end].min()
        area = df_t["dff"].iloc[start:end].sum()
        glucose = df_t["glucose"]
        event_glucose = glucose.iloc[start]
        #again this dataframe grew organically depending on what was needed this very second
        event_rows.append({
                "dataset": dataset,
                "condition": condition,
                "trace_id": tid,
                "event_id": event_id,
                "start_idx": start,
                "end_idx": end,
                "duration seconds": duration_s,
                "duration total": total_duration,
                "amplitude": amp,
                "event_area": area,
                "fs": fs,
                "glucose": event_glucose
        })
df_events = pd.DataFrame(event_rows)

#I am in fact aware of the fact that plotting and analysis should not be together in the same script
#I just also did not care at that exact moment
#In the "great restructuring event" that will happen at some point, hopefully, I will tidy this up
#but I really love these plots because they let me check right away where I messed up
#as for the rasterization and canvas.draw comments: this plot was needed for a poster
#and it refused to not be rasterized
#which messed up the scaling
#so I forced it
colors = ["red", "green", "blue", "orange"]
fig, ax = plt.subplots(figsize=(12, 6))
df_d1 = df_events[df_events["dataset"]=="dataset1"]
ax.set_rasterization_zorder(1)

for i, tid in enumerate(df_d1.trace_id.unique()):
   
    df_t = df_d1[df_d1.trace_id == tid]
    color = colors[i//4]
    for _, ev in df_t.iterrows():
        ax.plot([ev.start_idx, ev.end_idx], [i, i], color=color, linewidth=4, rasterized = True)

ax.set_yticks(range(len(df_d1.trace_id.unique())))
ax.set_yticklabels(df_d1.trace_id.unique())
ax.set_xlabel("Time index")
ax.set_ylabel("Trace")
ax.set_title("Event raster across all traces")
fig.canvas.draw()
fig.savefig("EventRaster_plateau2.svg")
plt.show()
plt.close()

#this one is more interesting because it tells me one specific onset
#I was playing with that idea for synchronization and it starts to look good
#therefore it gets to stay too
df_events["start_time"] = df_events["start_idx"]/df_events["fs"]
stim_start = 9800
stim_events = df_events[df_events["start_idx"] > stim_start]
onsets = stim_events.groupby("trace_id")["start_time"].min()
plt.hist(onsets, bins=20)
plt.xlabel("Onset time (s)")
plt.ylabel("Number of cells")
plt.title("Distribution of stimulus onset times")
plt.show()

print(onsets.describe())


#this is still doing debugging work because for some reason my Calcium peaks are being detected but not correctly saved
#therefore for voltage I get real numbers, for Calcium I get zero although in all later cases there are events
#this is interesting I still have to find a solution for it
for tid in df_traces.trace_id.unique():
    df_t = df_traces[df_traces.trace_id == tid]
    z = df_t["z"].values
    events = detect_events(z, fs, enter_threshold, exit_threshold, min_length, gap)
    print(tid, len(events))
#again hardcoded, hidden in the middle of the script, might overwrite things
#change at every iteration
#I will change it at some point I promise
df_events.to_pickle("df_events_plateau2.pkl")

summary_rows = []

for tid in df_traces.trace_id.unique():
     df_t = df_traces[df_traces.trace_id == tid]
     df_e = df_events[df_events.trace_id == tid]

     dataset = df_t["dataset"].iloc[0]
     condition = df_t["condition"].iloc[0]
     fs = df_t["fs"].iloc[0]
     glucose = df_t["glucose"]
     unique_glucose = np.unique(glucose)
     n_events = len(df_e)
     mean_amp = df_e["amplitude"].mean() if n_events > 0 else np.nan
     median_amp = df_e["amplitude"].median() if n_events > 0 else np.nan
     max_amp = df_e["amplitude"].max() if n_events > 0 else np.nan
     mean_duration = df_e["duration seconds"].mean() if n_events > 0 else np.nan
     median_duration = df_e["duration seconds"].median() if n_events > 0 else np.nan
     total_area = df_e["event_area"].sum() if n_events > 0 else 0
     time = df_t["time"].values
     total_time_s = df_t["time"].iloc[-1] - df_t["time"].iloc[0]
     event_freq = n_events/total_time_s if total_time_s > 0 else np.nan
     scales = np.arange(1, 128)
     sampling_dt = 0.003
     scale_min = 1
     scale_max = 128
     trace = df_t["dff"].values
     f_dom, dominance_ratio, z_score_f = dominant_frequency(trace, scales, sampling_dt, scale_min, scale_max)
     #print("tid:", tid, "scale_dom:", scale_dom2, "f_dom:", f_dom2)
     #print("new trace shape:", trace.shape)
     #print("new trace first 10:", trace[:10])
     #for the record I hate nested for loops but this one was necessary
     #at least if I want my glucose trace have any informational meaning for plots
     for g in unique_glucose:
        mask = (glucose == g)
        segment_signal = trace[mask]
        segment_time = time[mask]
        auc = np.trapz(segment_signal, segment_time)
        abs_auc = np.trapz(np.abs(segment_signal), segment_time)
        summary_rows.append({
            "dataset": dataset,
            "condition": condition,
            "trace_id": tid,
            "fs": fs,
            "n_events": n_events,
            "mean amplitude": mean_amp,
            "median amplitude": median_amp,
            "max amplitude": max_amp,
            "mean duration": mean_duration,
            "median duration": median_duration,
            "total area": total_area,
            "event frequency": event_freq,
            "dominant_frequency": f_dom,
            "dominance ratio": dominance_ratio,
            "frequency z-score": z_score_f,
            "glucose": g,
            "AUC": auc,
            "event Area": abs_auc
         })

df_summary = pd.DataFrame(summary_rows)
#at this point I am quite used to find my hidden hardcoded, data overwriting trip wires
#this one should be the last
#and it will be changed as well
df_summary.to_pickle("df_summary_plateau2.pkl")


#now these are basically just a lot of debugging plots
#and the point were I learnt about seaborn
#it is probably visible by the amount of seaborn suddenly used
plt.figure(figsize=(10, 6))

sns.swarmplot(
    data=df_summary,
    x="condition",
    y="dominant_frequency",
    size=6
)


plt.title("Dominant Frequency per Dataset")
plt.xlabel("condition")
plt.ylabel("Dominant Frequency (Hz)")
plt.show()

sns.boxplot(
    data = df_summary,
    x = "condition",
    y = "mean amplitude",
    width = 0.5
)
plt.show()

sns.boxplot(
    data = df_summary,
    x = "condition",
    y = "mean duration",
    width = 0.5
)
plt.show()

sns.boxplot(
    data = df_events,
    x = "dataset",
    y = "duration total",
    width = 0.5
)
plt.show()

plt.figure(figsize=(10, 6))

sns.barplot(
    data=df_summary,
    x="condition",
    y="total area",
    hue = "condition",
    palette=["#002D70", "#FF8000"]
)
#for anyone wondering #002D70 (pretty dark blue) and #FF8000 (orange) are just the colours of my university logo
#I tried to keep important plots in these colours for the poster to avoid too much happening at once
#whoever still remembers the raster plot also knows I failed at that
plt.show()

sns.swarmplot(
    data=df_summary,
    x="condition",
    y="total area",
    color="black",
    alpha=0.6
)

plt.title("Total Event Area per Dataset")
plt.xlabel("condition")
plt.ylabel("Total Event Area")
plt.show()

sns.swarmplot(data=df_summary, x="condition", y="event frequency", color = "black", alpha = 0.6)
plt.title("total event frequency per dataset")
plt.xlabel("condition")
plt.ylabel("Total Event Frequency")
plt.show()

#this figure also became important therefore the 600 dpi saving outcommented to not overwrite the actual nice figure
fig = plt.figure(figsize=(8, 4))
order = ["2mM", "10mM", "KCl"]
grouped = df_summary.groupby(["glucose", "condition"])["event Area"].mean().reset_index()
sns.barplot(
    data=grouped,
    x = "glucose",
    y = "event Area",
    hue = "condition",
    order = order,
    hue_order=["control", "ablated"],
    palette=["#002D70", "#FF8000"]
)
fig.canvas.draw()
#plt.savefig("barplot_glucose_platea.pdf", dpi = 600)
plt.show()
plt.show()

#this plot is the perfect debugging tool and horrible to run while running the script on a 8GB RAM laptop 
#running it on a better laptop makes it better
#anyways that is why it is there but outcommented it is basically the "peak detection debugging"-plot
'''
for dataset_name, df_d in df_traces.groupby("dataset"):
    trace_ids = df_d.trace_id.unique()
    n_cells = len(trace_ids)

    fig, axes = plt.subplots(
        n_cells, 1,
        figsize = (16, 7*n_cells),
        sharex = True
    )
    for ax, tid in zip(axes, trace_ids):
        df_t = df_traces[df_traces.trace_id == tid]
        z = np.array(df_t["z"])
        #z = df_t["z"]
        raw = df_t["signal"]
        dff_ex = df_t["dff"]

        df_e = df_events[df_events.trace_id == tid]
        merged_events = list(zip(df_e["start_idx"], df_e["end_idx"]))

        ax.plot(z, color="black")
        ax.axhline(enter_threshold, color="orange", linestyle = "--")
        ax.axhline(exit_threshold, color="green", linestyle = "--")
        for start, end in merged_events:
            ax.axvspan(start, end, color = "red", alpha = 0.2)
        ax.set_title(f"{dataset_name}-ROI-{tid}")
    plt.tight_layout()
    plt.show()
'''


