# Import Modules:
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
# Declare Variables:
assign_width = False
# Declare Functions:
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

# Declare classes:
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

    def find_n_derivative(self, n_derivatives):
        y = self.y
        n = n_derivatives
        all_ders = np.empty([n_derivatives, self.y.shape[0]])
        # der_y = np.empty_like(y)
        for i in range(n):
            y1 = np.gradient(y)
            all_ders[i] = y1
            y = y1
        self.derivatives = all_ders
    
    def plotting(self):
        fig = plt.figure()
        fig.canvas.set_window_title('Spectra Comparison.fig')
        plt.plot(self.x, self.y, color='red', linewidth=0.5, label='Raw Data')
        plt.plot(self.x, self.derivatives[0], color='blue', linewidth=0.5, label='First Derivative')
        plt.plot(self.x, self.derivatives[1], color='green', linewidth=0.5, label='Second Derivative')
        plt.legend()
        plt.show()


# Main Loop
pri_spectra = Spectra(file_import_selection())
pri_spectra.x_y_assignment()
pri_spectra.x_axis_limiting()
pri_spectra.find_n_derivative(2)
pri_spectra.plotting()

