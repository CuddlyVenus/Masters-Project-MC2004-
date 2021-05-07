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
assign_width = True # Boolean for assigning width
import_skip_rows = True # Boolean for skipping rows
normalise_data = True # Boolean for normalising y array
save_plot_figure = False # Boolean for saving the figure 

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
def coordinate_interpolateion(x, y):
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
            self.y = y1 # updates old y array with normalised y array 

    
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
        x, y = self.x, self.y
        # ---- Finding the first and second derivatives of the spectrum
        y_first_derivative = np.gradient(y)
        smoothed_y_first_derivative = signal.savgol_filter(y_first_derivative, 31, 2)
        y_second_derivative = np.gradient(smoothed_y_first_derivative)
        smoothed_y_second_derivative = signal.savgol_filter(y_second_derivative, 31, 2)
        self.y_first_derivative = smoothed_y_first_derivative * 1 
        self.y_second_derivative = smoothed_y_second_derivative * 100 # increasing y array for visability 
        # ---- First derivative peak detection
        zero_array = np.zeros_like(self.y)
        self.first_dev_intersection_indices = np.argwhere(np.diff(np.sign(zero_array - self.y_first_derivative))).flatten() # Plots zero line and determines all intercetions
        # against the y array
        for peak in self.first_dev_intersection_indices: # Loop to determine which zero intersections are have a postive gradient (maximas)
            lower_peak = peak - 1
            higher_peak = peak + 1
            lower_peak_y, peak_y, higher_peak_y = self.y_first_derivative[lower_peak], self.y_first_derivative[peak], self.y_first_derivative[higher_peak]
            arr = np.array([lower_peak_y, peak_y, higher_peak_y])
            bool_arr = (peak_y < lower_peak_y) & (peak_y > higher_peak_y) # Boolian response for postive gradient  
            if bool_arr == False:
                self.first_dev_intersection_indices = np.delete(self.first_dev_intersection_indices, np.where(self.first_dev_intersection_indices == peak)) # removes all 
                # negative gradient intersections form peak list

        # ---- Scipy peak detection for second derivation minima
        y2_second_derivative = self.y_second_derivative * -1 # inversion for peak minima detection 
        spectrum_peaks, _ = signal.find_peaks(y2_second_derivative, distance=10, prominence=0.01) # Values for peak detection are only best find for olivine atm
        # self.y_second_derivative = y_second_derivative * 100
        # ------ Scipy peak detection for second derivative maxima
        self.spectrum_peaks = spectrum_peaks
        self.all_peaks = np.concatenate((spectrum_peaks, self.first_dev_intersection_indices), axis=None)


    def peak_deconvolution(self):
        # ------ Rectangle based peak deconvolution
        # n = self.spectrum_peaks.size
        # self.all_peak_area_information = np.empty([n,4])
        # for i in range(n):
            # Generating the peak location at half height
            # peak = self.spectrum_peaks[i]
            # peak_height = self.y[peak]
            # half_peak_height = peak_height * 0.9 # set value for heigh of width determination 
            # Finding the intersecting points along the half height
            # hph_array = np.full_like(self.x, half_peak_height)
            # peak_intersection_indices = np.argwhere(np.diff(np.sign(hph_array - self.y))).flatten()
            # Sorting the intersecting points for the FWHM
            # peak_width_closest_intersection = self.x[peak_intersection_indices[(np.abs(peak_intersection_indices - peak)).argmin()]]
            # diff = np.abs(self.x[peak]-peak_width_closest_intersection)
            # peak_width_x_max = self.x[peak] + diff
            # peak_width_x_min = self.x[peak] - diff
            # peak_width_x_max = self.x[np.amin(peak_intersection_indices[peak_intersection_indices>peak])]
            # peak_width_x_min = self.x[np.amax(peak_intersection_indices[peak_intersection_indices<peak])]
            # Calculating peak area
            # peak_area_width = peak_width_x_max - peak_width_x_min
            # peak_area_height = half_peak_height
            # peak_area = peak_area_height * peak_area_width
            # Storing peak width and area information within a numpy array
            # self.all_peak_area_information[i] = np.array((peak_width_x_min,0,peak_area_width,peak_area_height))

        # ------ Limited Gaussian Peak Fitting 
        n = self.all_peaks.size
        i = self.y.size
        self.all_gussian_y_arrays = np.empty([n,i]) # generating empty array for deconvolved gaussian peaks 
        for i in range(n):
            peak = self.all_peaks[i]
            peak_height, peak_center = self.y[peak], self.x[peak]
            half_peak_height = peak_height * 0.8
            hph_array = np.full_like(self.x, half_peak_height) # intercetion y array
            peak_intersection_indices = np.argwhere(np.diff(np.sign(hph_array - self.y))).flatten()
            peak_width_closest_intersection = self.x[peak_intersection_indices[(np.abs(peak_intersection_indices - peak)).argmin()]]
            hwhm_sigma = (np.abs(self.x[peak]-peak_width_closest_intersection)) / 2.35482004503 # takes the FWHM and converts to sigma value (needs reference)
            # self.peak_y_gussian = peak_height*(1/(hwhm_sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((self.x-peak_center)/hwhm_sigma)**2))) # Reference Gussian calculation
            peak_y_gussian = peak_height * np.exp(-((self.x - peak_center)**2)/((2 * hwhm_sigma)**2)) # Generic gaussian curve 
            self.all_gussian_y_arrays[i] = peak_y_gussian

    def ir_band_determination(self):
        n = self.all_peaks.size
        for i in range(n):
            peak = self.all_peaks[i]
            peak_x, peak_y = self.x[peak], self.y[peak]
            # IR band specifications
            # ---- 1080 - 1175 wavenumber band
            v3_SiO4 = (peak_x > 1080) & (peak_x < 1175) # Boolian response 
            v1_SiO4 = (peak_x > 780) & (peak_x < 800)
            v2_SiO4 = (peak_x > 650) & (peak_x < 750)
            v4_SiO4 = (peak_x > 410) & (peak_x < 510)
    
    def peak_information(self):
        peak_indecies = self.all_peaks
        all_gussian_y_arrays = self.all_gussian_y_arrays
        x = self.x
        n = peak_indecies.size
        all_peak_areas = np.empty(n)
        all_peak_wavenumber = np.empty(n)
        for i in range(n):
            peak = peak_indecies[i]
            peak_y = all_gussian_y_arrays[i]
            peak_area = np.trapz(peak_y, x)
            all_peak_areas[i] = peak_area
            peak_wavenumber = self.x[peak]
            all_peak_areas[i] = peak_wavenumber
        all_peak_information = np.array([all_peak_areas, all_peak_wavenumber]).T
        output_df = pd.DataFrame(data=all_peak_information, columns=["Peak Wavenumber(cm-1)", "Peak Area"])
        output_df['Peak Y Values'] = all_gussian_y_arrays.tolist()
        output_df.to_csv(f'Saved Output\{self.file_index}_output.csv', sep='\t', index=False)
    
    def plot(self):
        # -- Whole spectrum figure
        fig = plt.figure(1)
        gs = fig.add_gridspec(3, hspace=0)
        axs = gs.subplots(sharex=True)
        axs[0].set_title(self.file_name, fontsize=15)
        # -- Specific Plots
        axs[0].plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data') # Plotting of the spectrum
        for y_gussian_array in self.all_gussian_y_arrays:
            axs[0].plot(self.x, y_gussian_array, linewidth=1, linestyle='dashdot', color='black') # Plotting fitted gussian curves
        # for peak in self.all_peaks:
            # x, y = self.x[peak], self.y[peak]
            # axs[0].text(x,y, '({})'.format(x))
        axs[1].plot(self.x, self.y_first_derivative, color='blue', linewidth=1.5, label='First Derivative') # plotting first dev on sep axs
        #axs[1].plot(self.x, np.zeros_like(self.x), color='black', linewidth=1) # visualising zero line intersections
        axs[2].plot(self.x, self.y_second_derivative, color='green', linewidth=1.5, label="Second Derivative")# Plotting the second derv on sep axs
        #axs[1].scatter(self.x[self.first_dev_intersection_indices], self.y_first_derivative[self.first_dev_intersection_indices], color='black') # Visual of second derv peak locations
        #axs[2].scatter(self.x[self.spectrum_peaks], self.y_second_derivative[self.spectrum_peaks], color='black') # Visual of second derv peak locations
        axs[0].scatter(self.x[self.spectrum_peaks], self.y[self.spectrum_peaks], color='green', marker=6, label='Second Derivative Peak') # Visual of second derv peak locations
        axs[0].scatter(self.x[self.first_dev_intersection_indices], self.y[self.first_dev_intersection_indices], color='blue', marker=7, label='First Derivative Peak') # Visual 
        # of first derv peak locations
        # -- Testing Plots
        # plt.hlines(self.hph, self.peak_width_x_min, self.peak_width_x_max, colors='black', linestyles='solid')
        # plt.scatter(self.x[self.peak_intersection_indices], self.y[self.peak_intersection_indices], color='blue')
        # for peak_area in self.all_peak_area_information:
            # rect = patches.Rectangle((peak_area[0], peak_area[1]), peak_area[2], peak_area[3], linewidth=1, linestyle='dashed', edgecolor='black', facecolor='None')
            # axs[0].add_patch(rect)
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
            fig.savefig(f'Saved Output/{self.file_index}', dpi=100)
        # -- Individual peak figure
        fig2, axs2 = plt.subplots(3,3)
        n = 0
        for a in range(3):
            for b in range(3):
                axs2[a, b].plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data')
                axs2[a, b].plot(self.x, self.all_gussian_y_arrays[n], linewidth=1, color='blue', label='Fitted Peak')
                axs2[a, b].set_xlim([np.max(self.x), np.min(self.x)])
                axs2[a, b].invert_xaxis()
                axs2[a, b].title.set_text(f'{np.round(self.x[self.all_peaks[n]], 0)} wavenumber(cm-1)')
                handles, labels = axs2[a, b].get_legend_handles_labels()
                n += 1
        fig2.legend(handles, labels, loc='upper center')
        peak_plotting_t.stop()
        plt.show()
# Main Loop - Processing the spectrum
print("---------- START ----------")
# TODO: Add method to init all functions within the class (automate)
pri_spectrum = Spectrum()
pri_spectrum.file_import_selection()
pri_spectrum.x_y_assignment()
pri_spectrum.normalise_data()
# pri_spectrum.baseline_correction()
pri_spectrum.x_axis_limiting()
pri_spectrum.peak_detection()
pri_spectrum.peak_deconvolution()
# pri_spectrum.ir_band_determination()
pri_spectrum.peak_information()
pri_spectrum.plot()
total_runtime_t.stop()
print("-----------END-----------")