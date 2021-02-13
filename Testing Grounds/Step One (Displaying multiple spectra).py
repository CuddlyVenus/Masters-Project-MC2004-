# Importing modules:
import os
import os.path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks 
# Assigment of vairables 
data_dir = "Project Data" # Directary of data to be used

# Functions:
def list_full_paths(directory): # List of files in project data with full path
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stiches directory to form path

def x_axis_limits(func_dataframe): # Detecting x axis max min for plotting axis scale
    max_axis = func_dataframe.max(axis=0)
    min_axis = func_dataframe.min(axis=0)
    x_max_rounded_value = round(max_axis[0], 0)
    x_min_rounded_value = round(min_axis[0], 0)
    # print(x_max_rounded_value, x_min_rounded_value) # Testing output
    return (x_max_rounded_value, x_min_rounded_value) 

def x_y_assigment(func_dataframe):
    x = func_dataframe.iloc[:, 0] # Assinging x value based on column index over name
    y = func_dataframe.iloc[:, 1] # Assigning y value based on column index over name
    # print("-----BREAK------")
    # clsprint(x,y)
    return(x, y)

def minima_peak_assigment(func_dataframe):
    y = func_dataframe.iloc[:, 1] # Assigning y value based on column index over name
    x = func_dataframe.iloc[:, 0] # Assinging x value based on column index over name
    y_neg = y*-1 # inverts the y values as spectra peaks face downwards
    minima = find_peaks(y_neg, height = -0.95, distance = 20) # Identifys minima as (x,y) co-ords 2D array
    min_pos = x[minima[0]] # 1d array of minima x values
    min_height = y_neg[minima[0]] # 1d array of minima y values
    min_height_neg = min_height*-1
    return(min_pos, min_height_neg)

# Main Loop
print("-----------------------------START----------------------------------")

print(list_full_paths(data_dir))
# Importing ASCII data formate into a workable panads dataframe
df = pd.read_csv(list_full_paths(data_dir)[0], sep='	', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df.columns = ["Wavenumber", "Transmitance"] # Naming Columns
df2 = pd.read_csv(list_full_paths(data_dir)[7], sep='	', header=None, usecols=[0, 1]) # Odd seperator, files without header and removes excess data
df2.columns = ["Wavenumber", "Transmitance"] # Naming Columns
# print(df.head(6)) # Test input

# Plotting Graph
fig = plt.figure()
fig.canvas.set_window_title('Spectra Comparison.fig')
(ax, ax1) = fig.subplots(2, sharex=True)
ax.set_title(list_full_paths(data_dir)[0])
ax1.set_title(list_full_paths(data_dir)[7])
ax.plot(x_y_assigment(df)[0], x_y_assigment(df)[1], color='red', linewidth=0.5) # Sub plot one of line x,y
ax1.plot(x_y_assigment(df2)[0], x_y_assigment(df2)[1], color='red', linewidth=0.5) # Sub plot two of line x2,y2
ax.scatter(minima_peak_assigment(df)[0], minima_peak_assigment(df)[1], color = 'black', s = 5, marker = 'X', label = 'Prominant Peaks') # Sub plot one minima identification
ax1.scatter(minima_peak_assigment(df2)[0], minima_peak_assigment(df2)[1], color = 'black', s = 5, marker = 'X') # Sub plot two minima identification
fig.legend()
fig.suptitle("FT-IR Spectra Comparison")
plt.xlabel('Wavenumber (cm-1)')
fig.text(0.04, 0.5, 'Transmitance', ha='center', va='center', rotation='vertical')
plt.xlim(x_axis_limits(df)) # defines the x axid limites by the datasets max and min values
plt.show()

print("------------------------------END-----------------------------------")