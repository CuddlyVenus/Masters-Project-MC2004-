# Importing modules:
import os
import os.path
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal, ndimage
# Assignment of vairables 
data_dir = "Project Data" # Directory of data to be used
gaussian = False

# Functions:
def list_full_paths(directory): # List of files in project data with full path
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stitches directory to form path

def x_axis_limits(func_dataframe): # Detecting x axis max min for plotting axis scale
    max_axis = func_dataframe.max(axis=0)
    min_axis = func_dataframe.min(axis=0)
    x_max_rounded_value = round(max_axis[0], 0)
    x_min_rounded_value = round(min_axis[0], 0)
    # print(x_max_rounded_value, x_min_rounded_value) # Testing output
    return (x_max_rounded_value, x_min_rounded_value) 

def x_y_assignment(func_dataframe):
    x = func_dataframe.iloc[:, 0] # Assinging x value based on column index over name
    y = func_dataframe.iloc[:, 1] # Assigning y value based on column index over name
    # print(x,y) Testing output
    return(x, y)

def minima_peak_assignment(func_dataframe):
    y = func_dataframe.iloc[:, 1] # Assigning y value based on column index over name
    x = func_dataframe.iloc[:, 0] # Assinging x value based on column index over name
    y_neg = y*-1 # inverts the y values as spectra peaks face downwards
    minima = signal.find_peaks(y_neg, height = -0.95, distance = 20) # Identifies minima as (x,y) co-ords 2D array
    min_pos = x[minima[0]] # 1d array of minima x values
    min_height = y_neg[minima[0]] # 1d array of minima y values
    min_height_neg = min_height*-1
    return(min_pos, min_height_neg)

def gaussian_smoothing(func_dataframe):
    if gaussian == True:
        x = ndimage.gaussian_filter1d(func_dataframe.iloc[:, 0], sigma=1) # smooths x values based on gaussian curve 
        y = ndimage.gaussian_filter1d(func_dataframe.iloc[:, 1], sigma=1) # smooths y values based on gaussian curve 
        # TODO: Need to determine cause of false peak identification when sigma = 1
        func_dataframe.insert(0, 'Gaussian Wavenumber', x, False) # updates dataframe with new gaussian values x axis (Cheap hack)
        func_dataframe.insert(1, 'Gaussian Transmittance', y, False)# updates dataframe with new gaussian values y axis (cheap hack)
        # print(x,y) Testing output
        return(x,y)
    else:
        return None 
        
# Main Loop
print("-----------------------------START----------------------------------")

print(list_full_paths(data_dir))
# Importing ASCII data formate into a workable panads dataframe
# TODO: Need to build reading and writing files to a Dataframe to a function that can be called (Automation)
df = pd.read_csv(list_full_paths(data_dir)[0], sep='	', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df.columns = ["Wavenumber", "Transmittance"] # Naming Columns
df2 = pd.read_csv(list_full_paths(data_dir)[7], sep='	', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df2.columns = ["Wavenumber", "Transmittance"] # Naming Columns
# print(df.head(5))  # Shows pre manipulation datframe as imported from files 
# print(df2.head(5)) # Shows pre manipulation datframe as imported from files 

# Manipulation of dataframe
gaussian_smoothing(df)
gaussian_smoothing(df2)

# Plotting Graph
fig = plt.figure()
fig.canvas.set_window_title('Spectra Comparison.fig')
(ax, ax1) = fig.subplots(2, sharex=True)
ax.set_title(list_full_paths(data_dir)[0])
ax1.set_title(list_full_paths(data_dir)[7])
ax.plot(x_y_assignment(df)[0], x_y_assignment(df)[1], color='red', linewidth=0.5) # Sub plot one of line x,y
ax1.plot(x_y_assignment(df2)[0], x_y_assignment(df2)[1], color='red', linewidth=0.5) # Sub plot two of line x2,y2
ax.scatter(minima_peak_assignment(df)[0], minima_peak_assignment(df)[1], color = 'black', s = 5, marker = 'X', label = 'Prominant Peaks') # Sub plot one minima identification
ax1.scatter(minima_peak_assignment(df2)[0], minima_peak_assignment(df2)[1], color = 'black', s = 5, marker = 'X') # Sub plot two minima identification
fig.legend()
fig.suptitle("FT-IR Spectra Comparison")
plt.xlabel('Wavenumber (cm-1)')
fig.text(0.04, 0.5, 'Transmittance', ha='center', va='center', rotation='vertical')
plt.xlim(x_axis_limits(df)) # defines the x axis limits by the datasets max and min values
plt.show()
# print(df.head(5))   # Shows post manipulation dataframe 
# print(df2.head(5)) # Shows post manipulation dataframe 

print("------------------------------END-----------------------------------")