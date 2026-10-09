import tifffile as tiff
import numpy as np
import cv2
from scipy.stats import linregress
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from tqdm import tqdm
import tkinter as tk
from tkinter import filedialog

def pick_folder():
    root = tk.Tk()
    root.withdraw()
    folder = filedialog.askdirectory(title = "Select folder to your tiffs")
    return folder
folder = pick_folder()
print("Selected folder: ", folder)

qc_root = os.path.join(folder, "QC_report")
qc_figures = os.path.join(qc_root, "QC_Figures")

tiff_files = [
    os.path.join(folder, f)
    for f in os.listdir(folder)
    if f.lower().endswith((".tif", ".tiff"))
]

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
os.makedirs(qc_root, exist_ok=True)
os.makedirs(qc_figures, exist_ok=True)
#here the Quality control calculation is carried out
def QC_data(signal_path, dataset_name, condition, frame_rate):
    with tiff.TiffFile(signal_path) as tif:
        pages = tif.pages

        motion_trace = []
        focus_trace = []
        expression_trace = []

        prev_frame = None

        for i, page in enumerate(pages):
            if i % 5 != 0:
                continue
            frame = page.asarray()
            frame_small = cv2.resize(frame, (0, 0), fx=0.25, fy = 0.25, interpolation=cv2.INTER_AREA)
            if prev_frame is not None:
                motion_trace.append(np.mean(np.abs(frame_small-prev_frame)))
            prev_frame = frame_small

            focus_trace.append(cv2.Laplacian(frame_small, cv2.CV_64F).var())

            expression_trace.append(frame_small.mean())
    motion_trace = np.array(motion_trace)
    focus_trace = np.array(focus_trace)
    expression_trace = np.array(expression_trace)
    #motion_trace = np.mean(np.abs(TimeSeries[1:]-TimeSeries[:-1]), axis=(1, 2))
    motion_score = motion_trace.mean()

    #focus_trace = [cv2.Laplacian(frame, cv2.CV_64F).var() for frame in TimeSeries]
    focus_score = np.mean(focus_trace)
    z_motion_trace = focus_trace/focus_trace.mean()
    #expression_trace = TimeSeries.mean(axis = (1, 2))
    expression_level = expression_trace.mean()

    bleaching_rate = linregress(np.arange(len(expression_trace)), expression_trace).slope
#normalized metrics for better comparison
    motion_norm = motion_score/expression_level
    focus_norm = focus_score/expression_level
    expression_norm = expression_level/np.max(expression_trace)
    bleaching_norm = (bleaching_rate/expression_level)*100
    z_motion_score = np.std(focus_trace)/focus_score
    
    return pd.DataFrame({
        "dataset": [dataset_name],
        "condition": [condition],
        "frame_rate": [frame_rate],
        "motion trace": [motion_trace],
        "focus trace": [focus_trace],
        "motion score": [motion_score],
        "focus_score": [focus_score],
        "expression trace": [expression_trace],
        "expression level": [expression_level],
        "bleaching": [bleaching_rate],
        "motion norm": [motion_norm],
        "focus norm": [focus_norm],
        "expression norm": [expression_norm],
        "bleaching norm": [bleaching_norm],
        "z_motion": [z_motion_score],
        "z_motion trace": [z_motion_trace]
    })
    
#datasets = [
#    ("1minute_voltron_cell_cameraStreaming.tif", "dataset1", "control", 100),
#    ("5ms_4x4binning_9min_islet4.tif", "dataset2", "control", 200),
#    #("180sec_4x4binning_5ms_fish6.tif", "dataset3", "control", 200),
#    #("180sec_4x4Binning_3ms_fish4_whathaveIDone_real.tif", "dataset4", "control", 300),
#    ("120sec_4x4binning_3ms_fish2.tif", "dataset5", "control", 300),
#    #("corrected_03_06_26_5minutes2mM_test.tif", "dataset6", "control", 100),
#    #("corrected_120sec_4x4binning_10ms_fish1.tif", "dataset7", "control", 100),
#    #("corrected_180sec_2x2binning_10ms_fish6.tif", "dataset8", "control", 100),
#    ("exvivo_f1.tif", "dataset9", "control", 100),
#    ("Islet1_2mM_30sec_10ms_trial3.tif", "dataset10", "control", 100)
#]

QC_list = []
for file in tqdm(tiff_files, desc = "Running QC"):
    QC_list.append(QC_data(
        signal_path = file,
        dataset_name = os.path.basename(file),
        condition = "control",
        frame_rate= 100
    ))
QC_metrics = pd.concat(QC_list, ignore_index = True)

QC_scores = QC_metrics.drop(columns=[
    "motion trace",
    "focus trace",
    "expression trace",
    "z_motion trace"
])

numeric_cols = [
    "motion score",
    "focus_score",
    "expression level",
    "bleaching",
    "motion norm",
    "focus norm",
    "expression norm",
    "bleaching norm",
    "z_motion"
]
QC_numeric = QC_scores[numeric_cols]

QC_threshold = {
    "motion score warning": 50,
    "motion score bad": 70,
    "focus score warning": 5,
    "focus score bad": 0.5,
    "expression warning": 0,
    "expression bad": 0,
    "bleaching warning": -0.1,
    "bleaching bad": -0.3,
    "z-motion warning": 0.8,
    "z-motion bad": 1
}
#this is very crude and also of course hardcoded
#I intent to change the hardcoded values with actual distribution values from my datasets
def label_dataset(row):
    if row["z_motion"] > QC_threshold["z-motion bad"]:
        return "bad"
    warning = False
    if row["motion norm"] > QC_threshold["motion score bad"]:
        return "bad"
    elif row["motion norm"] > QC_threshold["motion score warning"]:
        warning = True
    if row["focus norm"] < QC_threshold["focus score bad"]:
        return "bad"
    elif row["focus norm"] < QC_threshold["focus score warning"]:
        warning = True
    if row["expression norm"] < QC_threshold["expression bad"]:
        return "bad"
    elif row["expression norm"] < QC_threshold["expression warning"]:
        warning = True
    if row["bleaching norm"] < QC_threshold["bleaching bad"]:
        return "bad"
    elif row["bleaching norm"] < QC_threshold["bleaching warning"]:
        warning = True
    return "warning" if  warning else "good"

QC_scores["QC_label"] = QC_scores.apply(label_dataset, axis = 1)
#save the QC scores in a table to later compare several dataframes with each other
print(QC_scores)
save_csv = os.path.join(folder, "QC_report", "QC_Figures", "QC_metrics.csv")
QC_scores.to_csv(save_csv, index = False)
#row = QC_metrics.loc[1]   # pick dataset 0, or any index

for idx, row in QC_metrics.iterrows():
    fig, axes = plt.subplots(1, 4, figsize = (25, 4))
    fig.suptitle(f"QC_report overview - dataset {idx}")

    axes[0].plot(row["focus trace"])
    axes[0].set_title("Focus")
    axes[0].set_xlabel("Frame")
    axes[0].set_ylabel("Laplacian Variance")

    axes[1].plot(row["motion trace"])
    axes[1].set_title("Motion")
    axes[1].set_xlabel("Frame")
    axes[1].set_ylabel("Frame difference")

    axes[2].plot(row["expression trace"])
    axes[2].set_title("Bleaching")
    axes[2].set_xlabel("Frame")
    axes[2].set_ylabel("brightness")

    axes[3].plot(row["z_motion trace"])
    axes[3].set_title("Z-motion")
    axes[3].set_xlabel("Frame")
    axes[3].set_ylabel("blur")

    fig.tight_layout()
    save_path = os.path.join(folder, "QC_report", "QC_Figures", f"dataset_{idx}_QC.png")
    fig.savefig(save_path)
    plt.close(fig)


#check heatmap again normalization is weird
QC_norm = (QC_numeric-QC_numeric.mean(axis=0))/QC_numeric.std(axis=0)
save_heatmap = os.path.join(folder, "QC_report", "QC_Figures", "heatmap.png")
plt.figure(figsize = (12, 8))
sns.heatmap(QC_norm, annot=True, cmap="viridis", fmt=".2f",
    yticklabels=QC_scores["QC_label"])
plt.title("QC metrics heatmap")
plt.xlabel("QC metrics")
plt.ylabel("Dataset (QC label)")
plt.tight_layout()
plt.savefig(save_heatmap)
plt.show()