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
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from skimage import draw


def csv_to_mask(csv_path, image_shape):
    df = pd.read_csv(csv_path)

    roi_masks = []
    for shape_id, group in df.groupby("index", sort = False):
        group = group.sort_values("vertex-index")
        #axis highly dependent on how one draws their ROIs for Stardist it is 1 and 2 (on a 2D image), for polygon drawing on a 4d image it is 2 and 3
        #just check out the csv yourself
        coords = group[["axis-1", "axis-2"]].to_numpy() 
        mask = np.zeros(image_shape, dtype = bool)
        rr, cc = draw.polygon(
            coords[:, 0],
            coords[:, 1],
            shape = mask.shape
        )

        mask[rr, cc] = True
        roi_masks.append(mask)
    return roi_masks

def pick_tiff():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title = "Select tiff files",
        filetypes=[("TIFF files", "*.tif *.tiff")]
    )
    root.destroy()
    return file_path


def pick_csv():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title = "Select csv files",
        filetypes= [("CSV files", "*.csv")]
    )
    root.destroy()
    return file_path

tiff_path = pick_tiff()
print("Selected tiff: ", tiff_path)

csv_path = pick_csv()
print("Selected csv: ", csv_path)

df = pd.read_csv(csv_path)

print(df.head(10))
print(df.columns.tolist())

with tiff.TiffFile(tiff_path) as tif:
    TimeSeries = tif.asarray(out = 'memmap')
N = 20
preview = TimeSeries[:N].mean(axis=0)


#bg_masks = shapes_to_masks(bg_Layer)
signal_masks = csv_to_mask(csv_path, TimeSeries.shape[1:])

#print("Background ROIs: ", len(bg_masks))
print("signals ROIs: ", len(signal_masks))

#bg_traces = [TimeSeries[:, m].mean(axis=1) for m in bg_masks]
signal_traces = [TimeSeries[:, m].mean(axis=1) for m in signal_masks]

#np.save("Calcium_mnx_bg.npy", bg_traces)
np.save("Calcium_mnx_signal.npy", signal_traces)
#control whether the ROIs fit what I can see as a cell or not
fig, ax = plt.subplots(figsize = (10, 10))
ax.imshow(preview, cmap="gray")
for i, mask in enumerate(signal_masks):
    ax.contour(mask, levels = [0.5], colors = "red", linewidths = 1)
ax.set_title("overlay")
ax.axis("off")
plt.tight_layout()
plt.show()

