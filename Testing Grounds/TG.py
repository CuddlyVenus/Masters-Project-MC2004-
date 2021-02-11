import pandas as pd
import matplotlib.pyplot as plt

# importing HDF Testing 
# Imports from csv file, specifically with large space margin, with no headers and only using columns 1 and 0
df = pd.read_csv("Project Data\953.txt", sep='	', header=None, usecols=[0, 1])
df.columns = ["Wavenumber", "Transmitance"]
print(df.head(6))

# Set up axis determination functions
def max_x_axis(func_dataframe):
    max_axis = func_dataframe.max(axis=0)
    x_max_rounded_value = round(max_axis[0], 0)
    print(x_max_rounded_value)
    return x_max_rounded_value


def min_x_axis(func_dataframe):
    min_axis = func_dataframe.min(axis=0)
    x_min_rounded_value = round(min_axis[0], 0)
    print(x_min_rounded_value)
    return x_min_rounded_value


# Plotting Graph
ax1 = df.plot.line(x="Wavenumber", y="Transmitance", color='red', linewidth=0.5)
plt.xlim(min_x_axis(df), max_x_axis(df))
ax1.invert_xaxis()
ax1.plot()
plt.show()