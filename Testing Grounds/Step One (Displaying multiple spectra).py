# Importing modules:
import os
import os.path
import statistics
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, ndimage, sparse
from scipy.sparse.linalg import spsolve
from codetiming import Timer

# Assignment of code timer 
t = Timer(name="class", text="Total Runtime: {:.4f}s")

# Assignment of variables
t.start()
data_dir = "Project Data" # Directory of data to be used
smoothing = False # Selecting the type of smoothing to occur
bas_cor = False # Determining if a basline correction is to be applied to the spectrum
p = 0.001 # asymmetry assignment
lam = 102 # Smoothness assignment
Qat_dataframe_search = False # user spectra match or internal dataframe spectra search

# Functions:
def list_full_paths(directory): # List of files in project data with full path
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stitches directory to form path

def computeTicks (x, step = 100): # A function for assigning axis tick values based on a 1d array
    xMax, xMin = math.ceil(max(x)), math.floor(min(x))
    dMax, dMin = xMax + abs((xMax % step) - step) + (step if (xMax % step != 0) else 0), xMin - abs((xMin % step))
    return range(dMin, dMax, step)

def x_axis_limits(func_dataframe): # Detecting x axis max min for plotting axis scale
    max_axis = func_dataframe.max(axis=0)
    min_axis = func_dataframe.min(axis=0)
    x_max_rounded_value = round(max_axis[0], 0)
    x_min_rounded_value = round(min_axis[0], 0)
    return (x_max_rounded_value, x_min_rounded_value) 

def x_y_assignment(func_dataframe):
    x = np.array(func_dataframe.iloc[:, 0]) # Assigning x value based on column index over name
    y = np.array(func_dataframe.iloc[:, 1]) # Assigning y value based on column index over name
    return(x, y)

def maxima_peak_assignment(x_array, y_array):
    peaks, _ = signal.find_peaks(y_array, distance=20, prominence=0.01) 
    x_peaks = x_array[peaks]
    y_peaks = y_array[peaks]
    return(x_peaks, y_peaks)

def prominence_assignment(x_array, y_array):
    peaks, _ = signal.find_peaks(y_array, distance=20, prominence=0.01) 
    prominence, left_bases, right_bases = signal.peak_prominences(y_array, peaks)# produces the indices of the peaks in relation to y axis
    y_min = y_array[peaks] - prominence
    y_max = y_array[peaks]
    x = x_array[peaks]
    return(x, y_min, y_max)

def width_assignment(x_array, y_array):
    peaks, _ = signal.find_peaks(y_array, distance=20, prominence=0.01)
    prominence, left_bases, right_bases = signal.peak_prominences(y_array, peaks)# produces the indices of the peaks in relation to y axis
    # results_half = signal.peak_widths(y_array, peaks, rel_height=0.5, prominence_data=(prominence, left_bases, right_bases))
    results_full = signal.peak_widths(y_array, peaks, rel_height=1, prominence_data=((y_array[peaks] - prominence), left_bases, right_bases))
    y = y_array[peaks] - prominence
    x_min = results_full[2]
    x_max = results_full[3]
    return(y, x_array[left_bases], x_array[right_bases])

def graph_peak_arrangement(y_array): # Determins a mode baseline and assesses peak y direction, inverts for negative peak direction
    y_arr = np.array(y_array)
    mode_y_arr = statistics.mode(np.around(y_arr,decimals=2))
    higher_arr = 0
    lower_arr = 0
    for y in y_arr:
        if y >= mode_y_arr:
            higher_arr += 1
        if y < mode_y_arr:
            Lower_arr += 1
    if lower_arr > higher_arr:
        return(True)
    if higher_arr > lower_arr:
        return(False)
    else:
        return(None)

def baseline_correction(y, lam, p, niter=10): # "Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens in 2005
    if bas_cor == True:
        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
        w = np.ones(L)
        W = sparse.spdiags(w, 0, L, L)
        for i in range(niter):
            W.setdiag(w) # Do not create a new matrix, just update diagonal values
            Z = W + D
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return(z)
    else: 
        return(None)


def gaussian_smoothing(func_dataframe):
    if smoothing == True:
        x = ndimage.gaussian_filter1d(func_dataframe.iloc[:, 0], sigma=1) # smooths x values based on gaussian curve 
        y = ndimage.gaussian_filter1d(func_dataframe.iloc[:, 1], sigma=1) # smooths y values based on gaussian curve 
        # TODO: Need to determine cause of false peak identification when sigma = 1
        func_dataframe.insert(0, 'Gaussian Wavenumber', x, False) # updates dataframe with new gaussian values x axis (Cheap hack)
        func_dataframe.insert(1, 'Gaussian Transmittance', y, False)# updates dataframe with new gaussian values y axis (cheap hack)
        return(x,y)
    else:
        return None

def transpose_record_peaks(x_array, y_array): # Records and formates the peak coords within pandas dataframe
    a = pd.DataFrame(maxima_peak_assignment(x_array, y_array)) # Addition of x,y coord 1D array to panads dataframe
    b = a.transpose() # Transposes dataframe row to column for x,y titles
    return b # returns the final dataframe 
    #TODO: Find way to remove old row index values and replace with decending numericals, easy of access

def peak_comparison():
    #TODO: Compare data sets to determine if they are statistically similar or not and display result
    None

# Main Loop
print("-----------------------------START----------------------------------")

file_list = np.array(list_full_paths(data_dir))
# Importing ASCII data formate into a workable panads dataframe
# TODO: Need to build reading and writing files to a Dataframe to a function that can be called (Automation)
# TODO: Need to implement method for reading files with differing separators

# Primary dataframe designation and function assignment 
df = pd.read_csv(list_full_paths(data_dir)[0], sep=' ', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df.columns = ["Wavenumber", "Transmittance"] # Naming Columns
gaussian_smoothing(df)
pri_x,pri_y = x_y_assignment(df)
pri_x_peaks,pri_y_peaks = maxima_peak_assignment(pri_x, pri_y)
pri_pro_x, pri_pro_ymin, pri_pro_ymax = prominence_assignment(pri_x, pri_y)
pri_width_full_y, pri_width_full_x_min, pri_width_full_x_max = width_assignment(pri_x, pri_y)
pri_peakdf = transpose_record_peaks(pri_x, pri_y)
pri_blc = baseline_correction(pri_y, lam, p)

# Secondary dataframe designation and function assignment
df2 = pd.read_csv(list_full_paths(data_dir)[3], sep=' ', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df2.columns = ["Wavenumber", "Transmittance"] # Naming columns
gaussian_smoothing(df2)
sec_x,sec_y = x_y_assignment(df2)
sec_x_peaks,sec_y_peaks = maxima_peak_assignment(sec_x, sec_y)
sec_pro_x, sec_pro_ymin, sec_pro_ymax = prominence_assignment(sec_x, sec_y)
sec_width_full_y, sec_width_full_x_min, sec_width_full_x_max = width_assignment(sec_x, sec_y)
sec_peakdf2 = transpose_record_peaks(sec_x, sec_y)
sec_blc = baseline_correction(sec_y, lam, p)

# Plotting Graph
fig = plt.figure()
fig.canvas.set_window_title('Spectra Comparison.fig')
(ax, ax1) = fig.subplots(2, sharex=True, sharey=True)
plt.xticks(computeTicks(sec_x)) # sets the x ticks as defined by assigned function 
ax.set_title(list_full_paths(data_dir)[0])
ax1.set_title(list_full_paths(data_dir)[7])
ax.plot(pri_x, pri_y, color='red', linewidth=0.5) # Sub plot one of line x,y
ax1.plot(sec_x, sec_y, color='red', linewidth=0.5) # Sub plot two of line x2,y2
ax.scatter(pri_x_peaks, pri_y_peaks, color = 'black', s = 5, marker = 'X', label = 'Prominent Peaks') # Sub plot one minima identification
ax1.scatter(sec_x_peaks, sec_y_peaks, color = 'black', s = 5, marker = 'X') # Sub plot two minima identification
ax.vlines(pri_pro_x, pri_pro_ymin, pri_pro_ymax, color='black', linewidth=1, linestyle='dashed')
ax1.vlines(sec_pro_x, sec_pro_ymin, sec_pro_ymax, color='black', linewidth=1, linestyle='dashed')
ax.hlines(pri_width_full_y, pri_width_full_x_min, pri_width_full_x_max, color='green')
ax1.hlines(sec_width_full_y, sec_width_full_x_min, sec_width_full_x_max, color='green')
fig.legend()
fig.suptitle("FT-IR Spectra Comparison")
plt.xlabel('Wavenumber (cm-1)')
fig.text(0.04, 0.5, 'Reflectance', ha='center', va='center', rotation='vertical')
# plt.xlim(x_axis_limits(df)) # defines the x axis limits by the datasets max and min values
plt.ylim(0.00, 1.20)
# ax.invert_yaxis()
t.stop()
plt.show()
print("------------------------------END-----------------------------------")