# Import Modules
import os
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from scipy import signal, interpolate
from codetiming import Timer
# User Input Selection:
assign_width = False # Allows the user to limit the x axis array to specific minimum and maximum values (e.g 700 - 1200 cm-1): True/False
import_skip_rows = True # Lets the user designate the rows skipped to only access the spectrum data
normalise_data = True
save_plot_figure = False 

# Timer Designation - To test performance and determine areas of significant slowdown: 
peak_detection_t = Timer(name="class", text="Peak Detection Runtime: {:.4f}s")
peak_deconvolution_t = Timer(name="class", text="Peak Deconvolution Runtime: {:.4f}s")
peak_plotting_t = Timer(name="class", text="Plotting Runtime: {:.4f}s")
total_runtime_t = Timer(name="class", text="Total Runtime: {:.4f}s")
total_runtime_t.start()
peak_plotting_t.start()
# Function Designation:

# -- File Importing:
def list_full_paths(directory): # List of files in project data with relative full path (Project Data\\Spectrum.file)
    return [os.path.join(directory, file) for file in os.listdir(directory)] # takes file names from directory and stitches directory to form relative path

# -- Interpolation Function:
def linear_interpolateion(x, y):
    points = np.array([x, y]).T  # a (nbre_points x nbre_dim) array

    # Linear length along the line:
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)

    alpha = np.linspace(distance.min(), int(distance.max()), len(x))
    interpolator =  interpolate.interp1d(distance, points, kind='cubic', axis=0)
    interpolated_points = interpolator(alpha)

    out_x = interpolated_points.T[0]
    out_y = interpolated_points.T[1]

    return out_x, out_y

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
        # ---- Finding the second derivative of the spectrum
        # y1 = np.gradient(np.gradient(y))
        y_second_derivative = signal.savgol_filter(y, 31, 2, deriv=2)
        y_second_derivative = y_second_derivative * 100
        self.x_interpolated, self.y_second_derivative = linear_interpolateion(self.x, y_second_derivative)

        # ---- Scipy peak detection for second derivation minima
        y2_second_derivative = self.y_second_derivative * -1
        spectrum_peaks, _ = signal.find_peaks(y2_second_derivative, distance=10, prominence=0.01)
        # self.y_second_derivative = y_second_derivative * 100
        # ------ Scipy peak detection for second derivative maxima
        # Assigning spectra class peak indices
        self.spectrum_peaks = spectrum_peaks


    def peak_deconvolution(self):
        # ------ Rectangle based peak deconvolution
        n = self.spectrum_peaks.size
        self.all_peak_area_information = np.empty([n,4])
        for i in range(n):
            # Generating the peak location at half height
            peak = self.spectrum_peaks[i]
            peak_height = self.y[peak]
            half_peak_height = peak_height * 0.9 # set value for heigh of width determination 
            # Finding the intersecting points along the half height
            hph_array = np.full_like(self.x, half_peak_height)
            peak_intersection_indices = np.argwhere(np.diff(np.sign(hph_array - self.y))).flatten()
            # Sorting the intersecting points for the FWHM
            peak_width_closest_intersection = self.x[peak_intersection_indices[(np.abs(peak_intersection_indices - peak)).argmin()]]
            diff = np.abs(self.x[peak]-peak_width_closest_intersection)
            peak_width_x_max = self.x[peak] + diff
            peak_width_x_min = self.x[peak] - diff
            # peak_width_x_max = self.x[np.amin(peak_intersection_indices[peak_intersection_indices>peak])]
            # peak_width_x_min = self.x[np.amax(peak_intersection_indices[peak_intersection_indices<peak])]
            # Calculating peak area
            peak_area_width = peak_width_x_max - peak_width_x_min
            peak_area_height = half_peak_height
            peak_area = peak_area_height * peak_area_width
            # Storing peak width and area information within a numpy array
            self.all_peak_area_information[i] = np.array((peak_width_x_min,0,peak_area_width,peak_area_height))


        # ------ Limited Gaussian Peak Fitting (assumed FWHM)
        # peak_heights, peak_centers = self.y[self.peaks], self.x[self.peaks]
        # fwhm_sigma = 10
        # for peak_height, peak_center in zip(peak_heights, peak_centers):
        #     peak_y_gussian = peak_height*(1/(fwhm_sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((self.x-peak_center)/fwhm_sigma)**2)))
    
    def plot(self):
        # -- Whole spectrum figure
        fig = plt.figure(1)
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True)
        axs[0].set_title(self.file_name, fontsize=15)
        # -- Specific Plots
        axs[0].plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data') # Plotting of the spectrum
        axs[1].plot(self.x, self.y_second_derivative, color='green', linewidth=1.5, label="Second Derivative")# Plotting the second derv
        axs[1].scatter(self.x[self.spectrum_peaks], self.y_second_derivative[self.spectrum_peaks], color='black') # Visual of second derv peak locations
        axs[0].scatter(self.x[self.spectrum_peaks], self.y[self.spectrum_peaks], color='black') # Visual of second derv peak locations
        # -- Testing Plots
        # plt.hlines(self.hph, self.peak_width_x_min, self.peak_width_x_max, colors='black', linestyles='solid')
        # plt.scatter(self.x[self.peak_intersection_indices], self.y[self.peak_intersection_indices], color='blue')
        for peak_area in self.all_peak_area_information:
            rect = patches.Rectangle((peak_area[0], peak_area[1]), peak_area[2], peak_area[3], linewidth=1, linestyle='dashed', edgecolor='black', facecolor='None')
            axs[0].add_patch(rect)
        # -- Axis Naming
        fig.text(0.5, 0.04, 'Wavenumber (cm-1)', ha='center', va='center', fontsize=15)
        fig.text(0.06, 0.5, 'Reflectance', ha='center', va='center', rotation='vertical', fontsize=15)
        # -- Axis Size Limiting
        plt.xlim(np.max(self.x), np.min(self.x))
        # -- Axis Ticks Designation
        axs[0].xaxis.set_major_locator(ticker.MultipleLocator(100))
        for ax in axs:
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
        # -- Other
        fig.legend()
        if save_plot_figure == True:
            fig.savefig(f'Saved Plots/{self.file_index}', dpi=100)
        # -- Individual peak figure
        # fig2 = plt.figure(2)
        # gs = fig2.add_gridspec(3, 3)

        peak_plotting_t.stop()
        plt.show()
# Main Loop - Processing the spectrum
print("---------- START ----------")
pri_spectrum = Spectrum()
pri_spectrum.file_import_selection()
pri_spectrum.x_y_assignment()
pri_spectrum.normalise_data()
# pri_spectrum.baseline_correction()
pri_spectrum.x_axis_limiting()
pri_spectrum.peak_detection()
pri_spectrum.peak_deconvolution()
pri_spectrum.plot()
total_runtime_t.stop()
print("-----------END-----------")