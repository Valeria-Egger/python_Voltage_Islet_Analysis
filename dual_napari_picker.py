#dual napari ROI picker for simultanous Voltage and Calcium imaging

import napari
import skimage.draw as draw
from skimage.draw import polygon
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
from qtpy.QtCore import QTimer
#from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
#from qtpy.QtWidgets import QWidget, QVBoxLayout

root = tk.Tk()
root.withdraw()
def pick_tiff():
    return filedialog.askopenfilename(
        title = "Select tiff files",
        filetypes=[("TIFF files", "*.tif *.tiff")]
    )
print("Select Voltage tiff: ")
tiff_voltage = pick_tiff()
print("Select Calcium tiff: ")
tiff_Calcium = pick_tiff()
print("Selected tiff: ", tiff_voltage, tiff_Calcium)

with tiff.TiffFile(tiff_voltage) as tif:
    TimeSeries_voltage = tif.asarray(out = 'memmap')

with tiff.TiffFile(tiff_Calcium) as tif:
    TimeSeries_Calcium = tif.asarray(out = 'memmap')


#I choose a bit of a greater review because of artefacts in the system the dual pictures are taken with
N = 200
preview_voltage = TimeSeries_voltage[250:270].mean(axis=0)
preview_Calcium = TimeSeries_Calcium[250:270].mean(axis=0)

shape = preview_voltage.shape
viewer = napari.Viewer()

image_layer = viewer.add_image(
    preview_voltage,
    name = "voltage",
    colormap = "gray"
)

image_layer = viewer.add_image(
    preview_Calcium,
    name = "Calcium",
    colormap = "green",
    opacity = 0.4
)



bg_Layer = viewer.add_shapes(
    name = "background",
    shape_type = "polygon",
    edge_color = "yellow",
    face_color = "yellow",
    opacity = 0.3
)

voltage_Layer = viewer.add_shapes(
    name = "signals_voltage",
    shape_type = "polygon",
    edge_color = "red",
    face_color = "red",
    opacity = 0.3
)

Calcium_Layer = viewer.add_shapes(
    name = "signals_Calcium",
    shape_type = "polygon",
    edge_color = "green",
    face_color = "green",
    opacity = 0.3
)

plt.ion()
fig, ax = plt.subplots()


def poly_to_mask(poly, shape):
    rr, cc = polygon(poly[:, 1], poly[:, 0], shape)
    mask = np.zeros(shape, dtype = bool)
    mask[rr, cc] = True
    return mask
def extract_trace(TimeSeries, mask):
    return TimeSeries[:, mask].mean(axis=1)


def update_plot(event = None):
    ax.clear()

    if len(bg_Layer.data) > 0:
        bg_mask = poly_to_mask(bg_Layer.data[0], shape)
        bg_trace = extract_trace(TimeSeries_Calcium, bg_mask)
        ax.plot(bg_trace, label = "Background", color = "yellow")
    for i, poly in enumerate(voltage_Layer.data):
        mask_v = poly_to_mask(poly, shape)
        trace_v = extract_trace(TimeSeries_voltage, mask_v)
        ax.plot(trace_v, label = "voltage", color = "red")
    for i, poly in enumerate(Calcium_Layer.data):
        mask_c = poly_to_mask(poly, shape)
        trace_c = extract_trace(TimeSeries_Calcium, mask_c)
        ax.plot(trace_c, label = "Calcium", color = "green")
    ax.legend()
    fig.canvas.draw_idle()

#here I ran into problems with the plotting updating so much it froze my laptop so I had to make sure it updates once the user (me) is done
plot_timer = QTimer()
plot_timer.setSingleShot(True)
plot_timer.setInterval(1000)
plot_timer.timeout.connect(update_plot)

def request_update(event=None):
    plot_timer.start()

voltage_Layer.events.data.connect(request_update)
Calcium_Layer.events.data.connect(request_update)
bg_Layer.events.data.connect(request_update)

napari.run()
plt.ioff()


voltage_polys = voltage_Layer.data
Calcium_polys = Calcium_Layer.data
background_poly = bg_Layer.data

voltage_masks = []
Calcium_masks = []
bg_masks = []
voltage_traces = []
Calcium_traces = []
bg_traces = []

for i in range(len(voltage_polys)):
    mask_voltage = poly_to_mask(voltage_polys[i], shape)
    voltage_masks.append(mask_voltage)

    trace_voltage = extract_trace(TimeSeries_voltage, mask_voltage)
    voltage_traces.append(trace_voltage)

for i in range(len(Calcium_polys)):
    mask_Calcium = poly_to_mask(Calcium_polys[i], shape)
    Calcium_masks.append(mask_Calcium)

    trace_Calcium = extract_trace(TimeSeries_Calcium, mask_Calcium)
    Calcium_traces.append(trace_Calcium)

for i in range(len(background_poly)):
    mask_bg = poly_to_mask(background_poly[i], shape)
    bg_masks.append(mask_bg)

    trace_bg = extract_trace(TimeSeries_Calcium, mask_bg)
    bg_traces.append(trace_bg)


print("Voltage shape:", TimeSeries_voltage.shape)
print("Calcium shape:", TimeSeries_Calcium.shape)

with tiff.TiffFile(tiff_voltage) as tif:
    print("Voltage axes:", tif.series[0].axes)

with tiff.TiffFile(tiff_Calcium) as tif:
    print("Calcium axes:", tif.series[0].axes)

poly = Calcium_Layer.data[0]
mask = poly_to_mask(poly, TimeSeries_Calcium.shape[1:])

print("mask shape:", mask.shape)
print("image shape:", TimeSeries_Calcium.shape[1:])
print("ROI pixels:", mask.sum())

plt.figure()
plt.imshow(preview_Calcium, cmap="gray")
plt.imshow(mask, cmap="Reds", alpha=0.4)
plt.show()

calcium_trace = TimeSeries_Calcium[:, mask].mean(axis=1)

plt.figure()
plt.plot(calcium_trace)
plt.show()


print("Background ROIs: ", len(bg_masks))
print("voltage ROIs: ", len(voltage_masks))
print("Calcium ROIs: ", len(Calcium_masks))
#compatible with the rest of the analysis as the other ROI picker is too
np.save("islet1_2mM_bg.npy", bg_traces)
np.save("islet1_2mM_voltage.npy", voltage_traces)
np.save("islet1_2mM_Calcium.npy", Calcium_traces)