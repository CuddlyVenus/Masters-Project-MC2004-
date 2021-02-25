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
    """
    Computes domain with given step encompassing series x
    @ params
    x    - Required - A list-like object of integers or floats
    step - Optional - Tick frequency
    Here lies dark fuckery 
    """
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
    y = y_array # Assigning y value based on column index over name
    x = x_array # Assigning x value based on column index over name
    # y_neg = y*-1 # inverts the y values as spectra peaks face downwards
    minima = signal.find_peaks(y, height = 0.3, distance = 20) # Identifies minima as (x,y) co-ords 2D array
    min_pos = x[minima[0]] # 1d array of minima x values
    min_height = y[minima[0]] # 1d array of minima y values
    # min_height_neg = min_height*-1
    return(min_pos, min_height)

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

print(list_full_paths(data_dir))
# Importing ASCII data formate into a workable panads dataframe
# TODO: Need to build reading and writing files to a Dataframe to a function that can be called (Automation)
# TODO: Need to implement method for reading files with differing separators

# Primary dataframe designation and function assignment 
df = pd.read_csv(list_full_paths(data_dir)[0], sep=' ', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df.columns = ["Wavenumber", "Transmittance"] # Naming Columns
gaussian_smoothing(df)
pri_x,pri_y = x_y_assignment(df)
pri_min_pos,pri_min_height_neg = maxima_peak_assignment(pri_x, pri_y)
pri_peakdf = transpose_record_peaks(pri_x, pri_y)
pri_blc = baseline_correction(pri_y, lam, p)

# Secondary dataframe designation and function assignment
df2 = pd.read_csv(list_full_paths(data_dir)[3], sep=' ', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df2.columns = ["Wavenumber", "Transmittance"] # Naming columns
gaussian_smoothing(df2)
sec_x,sec_y = x_y_assignment(df2)
sec_min_pos,sec_min_height_neg = maxima_peak_assignment(sec_x, sec_y)
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
ax.scatter(pri_min_pos, pri_min_height_neg, color = 'black', s = 5, marker = 'X', label = 'Prominent Peaks') # Sub plot one minima identification
ax1.scatter(sec_min_pos, sec_min_height_neg, color = 'black', s = 5, marker = 'X') # Sub plot two minima identification
fig.legend()
fig.suptitle("FT-IR Spectra Comparison")
plt.xlabel('Wavenumber (cm-1)')
fig.text(0.04, 0.5, 'Reflectance', ha='center', va='center', rotation='vertical')
plt.xlim(x_axis_limits(df)) # defines the x axis limits by the datasets max and min values
plt.ylim(0.00, 1.20)
# ax.invert_yaxis()
t.stop()
plt.show()
print("------------------------------END-----------------------------------")