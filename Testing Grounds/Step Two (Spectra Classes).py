# ----- Modules ----- 
import os
import os.path
import statistics
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import rampy as rp
import peakutils
from scipy import signal, ndimage, sparse
from scipy.sparse.linalg import spsolve
from codetiming import Timer

# ----- Designate Timer -----
t = Timer(name="class", text="Total Runtime: {:.4f}s")
t.start()

# ----- Variable Assignment -----
assign_width = True
# ----- File Import -----
def list_full_paths(directory): # List of files in project data with full path
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stitches directory to form path

def file_import_selection():
    directory = "Project data"
    file_list = [os.path.join(directory, file) for file in os.listdir(directory)]
    file_index = np.nonzero(file_list)
    file_display = np.vstack((file_list, file_index)).T
    print(file_display)
    user_selected_file_index = int(input("Select file number for input:"))
    df = pd.read_csv(file_list[user_selected_file_index], delim_whitespace=True, header=None) # Odd seperator, files without header and removes excess data
    df.columns = ["Wavenumber", "Transmittance"] # Naming Columns
    return(df)

# ----- Class Structure -----
class Spectra():
    def __init__(self, df):
        self.df = df
    
    def x_y_assignment(self):
        x = np.array(self.df.iloc[:, 0]) # Assigning x value based on column index over name
        y = np.array(self.df.iloc[:, 1]) # Assigning y value based on column index over name
        self.x, self.y = x, y
    
    def x_axis_limiting(self):
        x = self.x
        y = self.y
        if assign_width == False:
            self.x = x
            self.y = y
        if assign_width == True:
            a, b  = input("Input desired range:").split()
            x_min, x_max = [int(x) for x in [a,b]]
            idx = np.where(np.logical_and(x>=x_min, x<=x_max))
            x1 = x[idx]
            y1 = y[idx]
            self.x = x1
            self.y = y1

    def peak_profile_generation(self):
        self.prominence, left_bases, right_bases = signal.peak_prominences(self.y, peaks)# produces the indices of the peaks in relation to y axis1
        y_min = self.y[peaks] - self.prominence
        y_max = self.y[peaks]
        x = self.x[peaks]
        min_pwidth = self.x[left_bases]
        max_pwidth = self.x[right_bases]
        peak_widths = max_width - min_width
        # Width cordinate generation
        # width = signal.peak_widths(self.y, peaks, rel_height=1, prominence_data=(y_min, left_bases, right_bases))
        self.peaks_x, self.peaks_y = peaks_x, peaks_y 
        self.prom_x, self.prom_y_min, self.prom_y_max, self.min_pwidth, self.max_pwidth = x, y_min, y_max, min_pwidth, max_pwidth

    def baseline_correction(self):
        self.baseline_values = peakutils.baseline(self.y)
        # self.y = self.y - self.baseline_values

    def plot(self):
        fig = plt.figure()
        fig.canvas.set_window_title('Spectra Comparison.fig')
        plt.plot(self.x, self.y, color='red', linewidth=1.5) # Sub plot one of line x,y
        # plt.plot(self.x, self.baseline_y, color='red', linewidth=1.5) # Sub plot one of line x,y
        plt.plot(self.x, self.baseline_values, color='blue', linewidth=0.5, label='Baseline Correction')
        # plt.plot(self.x, self.base_ploy, "-", color='grey', label="polynomial")
        plt.scatter(self.peaks_x, self.peaks_y, color = 'black', s = 5, marker = 'X', label = 'Prominent Peaks') # Sub plot one minima identification
        # plt.scatter(self.peaks_x, self.eval_height, color = 'black', s = 5, marker = 'X', label = 'Eval_height')
        # plt.scatter(self.intersect_x, self.intersect_y, color = 'orange', s = 5, marker = 'o', label = 'Intercetion')
        plt.vlines(self.prom_x, self.prom_y_min, self.prom_y_max, color='black', linewidth=1, linestyle='dashed', label='Peak Prominence')
        # plt.hlines(self.prom_y_min, self.min_pwidth, self.max_pwidth, color='green')
        fig.legend()
        fig.suptitle("FT-IR Spectra Comparison")
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel("Reflectance")
        plt.ylim(0.00, 1.20)
        plt.show()

# ----- Main Loop -----
print("---------- START ----------")
pri_spectra = Spectra(file_import_selection())
pri_spectra.x_y_assignment()
pri_spectra.x_axis_limiting()
pri_spectra.baseline_correction()
pri_spectra.peak_profile_generation()
pri_spectra.plot()
print("-----------END-----------")
# ----- End Timer -----
t.stop()