# Import Modules
import os
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import signal
from codetiming import Timer
# User Input Selection:
assign_width = False # Allows the user to limit the x axis array to specific minimum and maximum values (e.g 700 - 1200 cm-1): True/False
import_skip_rows = True # Lets the user designate the rows skipped to only access the spectrum data
normalise_data = True
save_plot_figure = True 

# Timer Designation - To test performance and determine areas of significant slowdown: 
peak_detection_t = Timer(name="class", text="Peak Detection Runtime: {:.4f}s")
peak_deconvolution_t = Timer(name="class", text="Peak Deconvolution Runtime: {:.4f}s")
peak_plotting_t = Timer(name="class", text="Plotting Runtime: {:.4f}s")
total_runtime_t = Timer(name="class", text="Total Runtime: {:.4f}s")
total_runtime_t.start()
# Function Designation:

# -- File Importing:
def list_full_paths(directory): # List of files in project data with relative full path (Project Data\\Spectrum.file)
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stitches directory to form relative path

# Class Designation
class Spectrum():
    
    def file_import_selection(self):
        directory = "Project data" # designating the project data folder 
        file_list = [os.path.join(directory, file) for file in os.listdir(directory)] # Array of all file paths
        file_index = np.nonzero(file_list) # Index of file list
        file_display = np.vstack((file_list, file_index)).T # Combined file list and index. Transformed vertically
        print(file_display)
        if import_skip_rows == False:
            user_selected_file_index = int(input("Select file number for input:")) # user selects index from file_display
            user_selected_file_name = file_list[user_selected_file_index]
            df = pd.read_csv(file_list[user_selected_file_index], delim_whitespace=True, header=None) # File is read into pandas dataframe based on user selected index 
            df.columns = ["Wavenumber", "Transmittance"] # Columns named
            self.df, self.file_index, self.file_name = df, user_selected_file_index, user_selected_file_name
        if import_skip_rows == True:
            user_selected_file_index = int(input("Select file number for input:")) # user selects index from file_display
            user_selected_file_name = file_list[user_selected_file_index]
            user_selected_skip_rows = int(input("Number of rows to be skipped:"))
            df = pd.read_csv(file_list[user_selected_file_index], delim_whitespace=True, header=None, skiprows=user_selected_skip_rows) # File is read into pandas dataframe based on user selected index 
            df.columns = ["Wavenumber", "Transmittance"] # Columns named
            self.df, self.file_index, self.file_name = df, user_selected_file_index, user_selected_file_name
    
    def x_y_assignment(self):
        x = np.array(self.df.iloc[:, 0]) # Assigning x value based on column index over name
        y = np.array(self.df.iloc[:, 1]) # Assigning y value based on column index over name
        self.x, self.y = x, y # Spectrum x, y values added to class data

    def normalise_data(self):
        y = self.y
        if normalise_data == False:
            self.y = y
        if normalise_data == True:
            y_max = np.max(y)
            y1 = y / y_max
            self.y = y1

    
    def baseline_correction(self):
        None
    
    def x_axis_limiting(self):
        x = self.x
        y = self.y
        if assign_width == False: # Whole spectrum is retainted 
            self.x = x
            self.y = y
        if assign_width == True: # User to designated range for spectrum plotting 
            a, b  = input("Input desired range:").split() # User to enter two values (Minimum and maximum)[Development needed for ease of use]
            x_min, x_max = [int(x) for x in [a,b]] 
            idx = np.where(np.logical_and(x>=x_min, x<=x_max)) # returns the indices within the range specified 
            x1 = x[idx]
            y1 = y[idx]
            self.x, self.y = x1, y1 # Class x,y values updated with user defined range
    
    def peak_detection(self):
        y = self.y
        # y1 = np.gradient(np.gradient(y))
        y_second_derivative = signal.savgol_filter(y, 31, 3, deriv=2)
        self.y_second_derivative = y_second_derivative * 100
    
    def peak_deconvolution(self):
        None
    
    def plot(self):
        # -- Specific Plots
        plt.plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data') # Plotting of the spectrum
        plt.plot(self.x, self.y_second_derivative, color='green', linewidth=1.5, label="Second Derivative")
        # -- Axis Naming
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel("Reflectance")
        # -- Axis Size Limiting
        plt.xlim(np.max(self.x), np.min(self.x))
        # -- Axis Ticks Designation
        ax = plt.axes()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
        # -- Other
        plt.title(self.file_name)
        plt.legend()
        if save_plot_figure == True:
            plt.savefig(f'Saved Plots/{self.file_index}', dpi=100)
        plt.show()
# Main Loop - Processing the spectrum
print("---------- START ----------")
pri_spectrum = Spectrum()
pri_spectrum.file_import_selection()
pri_spectrum.x_y_assignment()
pri_spectrum.normalise_data()
pri_spectrum.x_axis_limiting()
pri_spectrum.peak_detection()
pri_spectrum.plot()
total_runtime_t.stop()
print("-----------END-----------")