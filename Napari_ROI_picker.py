

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
import tkinter as tk
from tkinter import filedialog
#from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
#from qtpy.QtWidgets import QWidget, QVBoxLayout


def pick_tiff():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title = "Select tiff files",
        filetypes=[("TIFF files", "*.tif *.tiff")]
    )
    return file_path
tiff_path = pick_tiff()
print("Selected tiff: ", tiff_path)

#memory mapping to make the code executable on my machine, might not be needed at a laptop with more than 8GB RAM
with tiff.TiffFile(tiff_path) as tif:
    TimeSeries = tif.asarray(out = 'memmap')
#tiny preview, assumptions that data is already motion corrected
N = 20
preview = TimeSeries[:N].mean(axis=0)
#napari is used because napari gives me the opportunity to look at data the same way as FIJI
#also previous trials with automated ROI detection failed
#will revisit this idea later
viewer = napari.Viewer()
image_layer = viewer.add_image(preview)

#plot_widget = QWidget()
#plot_layout = QVBoxLayout(plot_widget)
#fig, ax = plt.subplots()
#canvas = FigureCanvasQTAgg(fig)

#plot_layout.addWidget(canvas)
#viewer.window.add_dock_widget(plot_widget, area = "right")



bg_Layer = viewer.add_shapes(
    name = "background",
    shape_type = "polygon",
    edge_color = "yellow",
    face_color = "yellow",
    opacity = 0.3
)

signal_Layer = viewer.add_shapes(
    name = "signals",
    shape_type = "polygon",
    edge_color = "red",
    face_color = "red",
    opacity = 0.3
)
#interactive Matplotlib-plot to see how the traces look like, scan for interesting stuff to appear
plt.ion()
fig, ax = plt.subplots()

#within this plot define quiet and signal windows, for calculations later on
quiet_window = {}
signal_window = {}

def onselect_quiet(qmin, qmax):
    quiet_window["qmin"] = int(qmin)
    quiet_window["qmax"] = int(qmax)
    print("Quiet window: ", quiet_window)

def onselect_signal(smin, smax):
    signal_window["smin"] = int(smin)
    signal_window["smax"] = int(smax)
    print("signal Window: ", signal_window)

span_quiet = SpanSelector(
    ax,
    onselect_quiet,
    "horizontal",
    useblit = True,
    props=dict(alpha=0.3, facecolor="blue"),
    button = 1)

span_signal = SpanSelector(
    ax,
    onselect_signal,
    "horizontal",
    useblit = True,
    props=dict(alpha=0.3, facecolor="green"),
    button = 3)

lines = {}


#here we are back at the napari code
def compute_mask(shape):
    coords = np.array(shape)
    mask = np.zeros(TimeSeries.shape[1:], dtype = bool)
    rr, cc = draw.polygon(coords[:, 0], coords[:, 1], mask.shape)
    mask[rr, cc] = True
    return mask

#this actually does do the interactive Figure
def update_plot(event=None):
    ax.clear()

    if len(bg_Layer.data) > 0:
        bg_mask = compute_mask(bg_Layer.data[0])
        bg_Trace = TimeSeries[:, bg_mask].mean(axis=1)
        ax.plot(bg_Trace, label = "Background", color = "yellow")

    for i, shape in enumerate(signal_Layer.data):
        sig_mask = compute_mask(shape)
        sig_Trace = TimeSeries[:, sig_mask].mean(axis=1)
        ax.plot(sig_Trace, label = f"signal{i}", alpha = 0.8)

    ax.legend()
    fig.canvas.draw_idle()

#here the quiet and signal windows are important for calculations
def compute_snr():
    if "qmin" not in quiet_window or "smin" not in signal_window:
        print("Select a quiet and a signal window")
        return

    sig_mask = compute_mask(signal_Layer.data[0])
    sig_Trace = TimeSeries[:, sig_mask].mean(axis = 1)

    bg_mask = compute_mask(bg_Layer.data[0])
    bg_Trace = TimeSeries[:, bg_mask].mean(axis=1)

    qmin, qmax = quiet_window["qmin"], quiet_window["qmax"]
    smin, smax = signal_window["smin"], signal_window["smax"]

    quiet = sig_Trace[qmin: qmax]
    signal = sig_Trace[smin: smax]
    background = bg_Trace[:]
    signal_for_bg = sig_Trace[:]

    noise = np.std(quiet)
    baseline = np.mean(quiet)
    amplitude = np.min(signal)-baseline
    snr = amplitude/noise

    print("noise: ", noise)
    print("baseline: ", baseline)
    print("amplitude: ", amplitude)
    print("snr: ", snr)
    print("signal min/max:", signal_for_bg.min(), signal_for_bg.max())
    print("background min/max:", background.min(), background.max())
    

bg_Layer.events.data.connect(update_plot)
signal_Layer.events.data.connect(update_plot)

@magicgui(call_button="Compute SNR")
def snr_button():
    compute_snr()

viewer.window.add_dock_widget(snr_button, area="right")

napari.run()
plt.ioff()

#then actually extract the traces
def shapes_to_masks(shapes_layer):
    roi_masks = []
    for shape in shapes_layer.data:
        coords = np.array(shape)
        mask = np.zeros(TimeSeries.shape[1:], dtype = bool)
        rr, cc = draw.polygon(coords[:, 0], coords[:, 1], mask.shape)
        mask[rr, cc] = True
        roi_masks.append(mask)
    return roi_masks
    
bg_masks = shapes_to_masks(bg_Layer)
signal_masks = shapes_to_masks(signal_Layer)

print("Background ROIs: ", len(bg_masks))
print("signals ROIs: ", len(signal_masks))

bg_traces = [TimeSeries[:, m].mean(axis=1) for m in bg_masks]
signal_traces = [TimeSeries[:, m].mean(axis=1) for m in signal_masks]
#hardcoded, change for each new dataset
np.save("voltage1_bg.npy", bg_traces)
np.save("voltage1_signal.npy", signal_traces)