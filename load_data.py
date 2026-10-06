import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
import pylab

#actual script do not touch until you actually want to change something
#no exploratory science in here

#thresholds are still hardcoded because they are highly variable
#I want to revisit that decision later
enter_threshold = -3.0
exit_threshold = -2

def load_dataset(signal_path, background_path, dataset_name, condition, frame_rate, perfusion):
    signal_traces = np.load(signal_path)
    print(signal_traces.shape)
    background_traces = np.load(background_path)
    background_traces = np.squeeze(background_traces)

    n_traces, n_timepoints = signal_traces.shape
    time = np.arange(n_timepoints)/frame_rate
    #also every single perfusion version is hardcoded
    #I do not plan to visit that decision later    

    #hypoglycemic perfusion
    if perfusion == 1:
        glucose_trace = np.array([
            "2mM" if t < 180 else
            "1mM" if t < 360 else
            "0.5mM" if t < 540 else
            "0.25mM" if t < 720 else
            "0mM" if t < 900 else
            "KCl" 
            for t in range(len(time))
        ])
        #hyperglycemic perfusion
    elif perfusion == 2:
        glucose_trace = np.array([
                "2mM" if t < 180 else
                "6mM" if t < 360 else
                "9mM" if t < 540 else
                "12mM" if t < 720 else
                "15mM" if t < 900 else
                "18mM" if t < 1080 else
                "KCl" 
                for t in range(len(time))
            ])
        #most used voltage perfusion
    elif perfusion == 3:
        glucose_trace = np.array([
            "2mM" if t < 3000 else
            "10mM" if t < 9000 else
            "KCl"
            for t in range(len(time))
        ])
        #one of my experiments perfusion
    elif perfusion == 4:
        glucose_trace = np.array([
            "2mM" if t < 300 else
            "5mM"
            for t in range(len(time))
        ])
        #no perfusion, normoglycemia
    else:
        glucose_trace = np.array([
            "2mM"
            for t in range(len(time))
        ])
    rows = []

    for trace_index in range(n_traces):
        trace = signal_traces[trace_index]
        #actually this one is due to the fact that halfway through I realized most of the script can be reused for Calcium imaging too
        #so I had to distinguish the two with something
        #which in this case is framerate
        #everything above 100 is Voltage, everything below is Calcium
        #I will totally revisit this decision later because no
        if frame_rate > 99:
            cutoff = 0.05 #depends on Hz (300Hz 0.01, 100Hz 0.05)
        
            order = 2
            sos = butter(order, cutoff/(frame_rate/2), btype='low', output='sos')
            #sos_high = butter(order, 50/(frame_rate/2), btype='high', output='sos')
            baseline_trace = sosfiltfilt(sos, signal_traces[trace_index])
            #noise_trace = sosfiltfilt(sos_high, signal_traces[trace_index])
            dff_trace = (signal_traces[trace_index]-baseline_trace)/baseline_trace

            residual = dff_trace+np.median(dff_trace)
            quiet_half = residual[residual < np.percentile(residual, 50)]
            mad = np.median(np.abs(quiet_half))
            noise_std = 1.4826*mad
            #noise_std = np.clip(noise_std, 0.02, 0.2)
            z_trace = (dff_trace-np.median(dff_trace))/noise_std
        else:
            #peaceful min max scaling for Calcium imaging because they can deal with that
            baseline_trace = np.min(trace)
            dff_trace = (trace-baseline_trace)/(np.max(trace)-baseline_trace)
            z_trace = (dff_trace-np.mean(dff_trace))/np.std(dff_trace)
        #this dataframe grew organically depending on what I thought I needed right this second
        rows.append(pd.DataFrame({
         "dataset": dataset_name,
         "trace_id": f"{dataset_name}_{trace_index}",
         "condition": condition,
         "fs": frame_rate,
         "time": time,
         "signal": signal_traces[trace_index],
         "background": background_traces[0],
         "baseline": baseline_trace,
         "dff": dff_trace,
         "z": z_trace,
         "glucose": glucose_trace
        }))


    return pd.concat(rows, ignore_index=True)

#here I can make a list for all the data I have that I want to put in the dataframe
#here is just everything hardcoded
#I am planning on restructuring that as well at some point
datasets = [
    ("voltage_plateau2_signal.npy", "voltage_plateau2_bg.npy", "dataset1", "control", 100, 3),
    ("voltageposter2_signal.npy", "voltageposter2_bg.npy", "dataset2", "control", 100, 3), 
]



df_traces = pd.concat([
    load_dataset(*d) for d in datasets
], ignore_index = True)
#this saving line is not only hardcoded but also hidden in the middle of the script
#if you do not want to overwrite anything change at every iteration
#with some luck I will revisit that decision later on as well
df_traces.to_pickle("df_traces_plateau2.pkl")

#2 plots to check whether what we intended to save is actually being saved
for tid in df_traces.trace_id.unique():
    df_t = df_traces[df_traces.trace_id == tid]
    plt.plot(df_t["time"], df_t["z"], label=f"trace {tid}")

plt.legend()
plt.show()


df1 = df_traces[df_traces["dataset"]== "dataset2"]
trace_ids = df1.trace_id.unique()
n = len(trace_ids)

fig, axes = plt.subplots(n, 1, figsize = (10, 3*n), sharex = True)
for ax, tid in zip(axes, trace_ids):
    df_t = df1[df1.trace_id == tid]
    ax.plot(df_t["time"], df_t["z"])
    ax.set_title(f"Trace {tid}")
plt.tight_layout()
plt.show()
