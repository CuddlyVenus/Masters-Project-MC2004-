# Importing modules:
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks 
# Terminal break
print("-----------------------------START----------------------------------")
# Importing ASCII data formate into a workable panads dataframe
df = pd.read_csv("Project Data\954.txt", sep='	', header=None, usecols=[0, 1]) # Uniqe seperator, files without header and removes excess data
df.columns = ["Wavenumber", "Transmitance"] # Naming Columns
df2 = pd.read_csv("Project Data\960.txt", sep='	', header=None, usecols=[0, 1]) # Uniqe seperator, files without header and removes excess data
df2.columns = ["Wavenumber", "Transmitance"] # Naming Columns
# print(df.head(6)) # Test input

# Functions: 
def x_axis_limits(func_dataframe): # Detecting x axis max min for plotting axis scale
    max_axis = func_dataframe.max(axis=0)
    min_axis = func_dataframe.min(axis=0)
    x_max_rounded_value = round(max_axis[0], 0)
    x_min_rounded_value = round(min_axis[0], 0)
    # print(x_max_rounded_value, x_min_rounded_value) # Testing output
    return (x_max_rounded_value, x_min_rounded_value) 

# Defining the x and y arrys for plotting line spectra
x = df["Wavenumber"]
y = df["Transmitance"]
inversey = y*-1 # creating an inverse of the y data for peak identification (opposite to peaks on spectra)
x2 = df2["Wavenumber"]
y2 = df2['Transmitance']
inversey2 = y2*-1 # creating an inverse of the y data for peak identification (opposite to peaks on spectra)

# Finding the minima with relation to the y array
minima = find_peaks(inversey, height = -0.95, distance = 20) # Identifys minima as (x,y) co-ords 2D array
min_pos = x[minima[0]] # 1d array of minima x values
min_height = inversey[minima[0]] # 1d array of minima y values

minima2 = find_peaks(inversey2, height = -0.95, distance = 20)
min_pos2 = x2[minima2[0]]
min_height2 = inversey2[minima2[0]]


# Plotting Graph
fig = plt.figure()
fig.canvas.set_window_title('Instert File Name Here')
(ax, ax1) = fig.subplots(2, sharex=True)
ax.plot(x,y, color='red', linewidth=0.5) # Sub plot one of line x,y
ax1.plot(x2,y2, color='red', linewidth=0.5) # Sub plot two of line x2,y2
ax.scatter(min_pos, min_height*-1, color = 'black', s = 5, marker = 'X', label = 'Prominant Peaks') # Sub plot one minima identification
ax1.scatter(min_pos2, min_height2*-1, color = 'black', s = 5, marker = 'X') # Sub plot two minima identification
fig.legend()
fig.suptitle("Comparison")
plt.xlim(x_axis_limits(df)) # defines the x axid limites by the datasets max and min values
plt.show()