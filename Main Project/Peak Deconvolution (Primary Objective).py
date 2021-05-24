# Import Modules
import os
import os.path
import math
import numpy as np
from numpy.core.numeric import cross
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
peak_deconvolution_type = 'second'

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
# -- Curve Result Class:
class Curve_Results:

    def __init__(self, curve):
        self.curve = curve
        self.internal_comp = {}
        self.internal_comp_highest_value = []

    def FindHighestIntersectionArea(self):
        key_array = [*self.internal_comp]
        highest_value_array = []
        for key in self.internal_comp:
            max_value = np.max(self.internal_comp[key])
            highest_value_array.append(max_value)
        highest_value_array2 = np.array(highest_value_array)
        HV_index = np.argmax(highest_value_array2)
        self.internal_comp_highest_value.append(key_array[HV_index])
        self.internal_comp_highest_value.append(highest_value_array2[HV_index])

# -- Curve Class:
class Curve:

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def GetTotalArea(self):
        return math.sqrt(math.pi) * math.sqrt(2) * (self.a * self.c)

    def GetDefiniteArea(self, xmin, xmax):
        p1 = math.erf((xmax - self.b) / (math.sqrt(2) * self.c))
        p2 = math.erf((xmin - self.b) / (math.sqrt(2) * self.c))
        return (math.sqrt(math.pi) / math.sqrt(2)) * self.a * self.c * (p1 - p2)

    def GetIntersectionArea(self, other):
        intersections = self.GetIntersectionPoints(other)
        if (len(intersections) == 2):
            double_intersection = self
            xmin = min(intersections)
            xmax = max(intersections)
            area = double_intersection.GetDefiniteArea(xmin, xmax)
            return area

        elif (len(intersections) == 1):
            left = self if self.b <= other.b else other
            right = self if self.b > other.b else other
            left_area = left.GetDefiniteArea(intersections[0], 10000)
            right_area = right.GetDefiniteArea(-10000, intersections[0])
            return left_area + right_area

        else:
            peak_1 = self
            one_half_sigma1 = peak_1.c * 2.5
            peak1_sigma_range_min, peak1_sigma_range_max = (peak_1.b - one_half_sigma1), (peak_1.b + one_half_sigma1)
            peak_2 = other
            one_half_sigma2 = peak_2.c * 2.5
            peak2_sigma_range_min, peak2_sigma_range_max = (peak_2.b - one_half_sigma2), (peak_2.b + one_half_sigma2)
            if peak1_sigma_range_min > peak2_sigma_range_min and peak1_sigma_range_max < peak2_sigma_range_max:
                lesser = peak_1
                lesser_area = lesser.GetTotalArea()
                return lesser_area
            elif peak2_sigma_range_min > peak1_sigma_range_min and peak2_sigma_range_max < peak1_sigma_range_max:
                lesser = peak_2
                lesser_area = lesser.GetTotalArea()
                return lesser_area
            else:
                return 0

    def GetIntersectionPointsOLD(self, other):
        x = 1 / (2 * self.c ** 2) - 1 / (2 * other.c ** 2)
        y = other.b / (other.c**2) - self.b/(self.c**2)
        z = self.b**2 /(2*self.c**2) - other.b**2 / (2*other.c**2) - np.log((other.c*self.a)/(self.c*other.a))
        return np.roots([x,y,z])

    def GetIntersectionPoints(self, other):
        # Old method (above) is based on a normal PDF where a is proportional to c, so it doesn't work for a few cases.
        # below is based on the equation we are using, worked out using an online root finder. I think it works.
        a1 = self.a
        a2 = other.a
        b1 = self.b
        b2 = other.b
        c1 = self.c
        c2 = other.c
        if (c1 == c2):
            c2 += 0.000001
        # Change this when in report
        cum = (2 * c1**2) * math.log(a2) - (2 * c1 ** 2) * math.log(a1) + (b1 ** 2) + (b2 ** 2) + (2 * c2 ** 2) * math.log(a1) - 2 * b1 * b2 - (2 * c2**2) * math.log(a2)
        if cum < 0:
            return []

        common = c1 * c2 * math.sqrt(cum)
        x1 = (-(b1 * c2 ** 2) + (b2 * c1 ** 2) + common) / (c1 ** 2 - c2 ** 2)
        x2 = (-(b1 * c2 ** 2) + (b2 * c1 ** 2) - common) / (c1 ** 2 - c2 ** 2)

        ret = [x1, x2]
        return [x for x in ret if self.GetY(x) > self.a * 0.001 and other.GetY(x) > other.a * 0.001]

    def GetY(self, x):
        return self.a * np.exp(-((x - self.b) ** 2)/(2 * (self.c ** 2)))
# -- Spectrum Class:
class Spectrum():

    def __init__(self, name):
        self.name = name
    
    def file_import_selection(self):
        directory = "Project data" # designating the project data folder 
        file_list = [os.path.join(directory, file) for file in os.listdir(directory)] # Array of all file paths
        file_index = np.nonzero(file_list) # Index of file list
        file_display = np.vstack((file_list, file_index)).T # Combined file list and index. Transformed vertically
        print(file_display)
        if import_skip_rows == False:
            user_selected_file_index = int(input(f"Select file number for {self.name} spectrum:")) # user selects index from file_display
            user_selected_file_name = file_list[user_selected_file_index]
            df = pd.read_csv(file_list[user_selected_file_index], delim_whitespace=True, header=None) # File is read into pandas dataframe based on user selected index 
            df.columns = ["Wavenumber", "Transmittance"] # Columns named
            user_selected_file_name = user_selected_file_name.split('\\')[-1]
            self.df, self.file_index, self.file_name = df, user_selected_file_index, user_selected_file_name
        if import_skip_rows == True:
            user_selected_file_index = int(input(f"Select file number for {self.name} spectrum:")) # user selects index from file_display
            user_selected_file_name = file_list[user_selected_file_index]
            user_selected_skip_rows = int(input("Number of rows to be skipped:"))
            df = pd.read_csv(file_list[user_selected_file_index], delim_whitespace=True, header=None, skiprows=user_selected_skip_rows) # File is read into pandas dataframe based on user selected index 
            df.columns = ["Wavenumber", "Transmittance"] # Columns named
            user_selected_file_name = user_selected_file_name.split('\\')[-1]
            self.df, self.file_index, self.file_name = df, user_selected_file_index, user_selected_file_name
    
    def x_y_assignment(self):
        #TODO Need to force all x-axis to go from smallest to largeest value after import
        x = np.array(self.df.iloc[:, 0]) # Assigning x value based on column index over name
        y = np.array(self.df.iloc[:, 1]) # Assigning y value based on column index over name
        if x[0] > x[1]:
            x = np.flip(x)
            y = np.flip(y)
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
            a, b  = input(f"Input desired x range for {self.name}:").split() # User to enter two values (Minimum and maximum)[Development needed for ease of use]
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
        self.fist_dev_peaks = self.first_dev_intersection_indices
        if peak_deconvolution_type == 'first':
            self.all_peaks = self.fist_dev_peaks
        if peak_deconvolution_type == 'second':
            self.all_peaks = self.spectrum_peaks
        if peak_deconvolution_type == 'both':
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
        self.all_peak_areas = np.empty(n)
        self.all_peak_wavenumbers = np.empty(n)
        self.all_gussian_components = np.empty([n,3], dtype=np.object)
        for i in range(n):
            peak = self.all_peaks[i]
            peak_height_a, peak_center_b = self.y[peak], self.x[peak]
            half_peak_height = peak_height_a * 0.9
            hph_array = np.full_like(self.x, half_peak_height) # intercetion y array
            # if peak intersection indices is empty then skip:
            peak_intersection_indices = np.argwhere(np.diff(np.sign(hph_array - self.y))).flatten()
            if peak_intersection_indices.size == 0:
                peak_intersection_indices = [0, (self.y.size-1)]
            peak_width_closest_intersection = self.x[peak_intersection_indices[(np.abs(peak_intersection_indices - peak)).argmin()]]
            hwhm_sigma_c = (np.abs(self.x[peak]-peak_width_closest_intersection)) / 2.35482004503 # takes the FWHM and converts to sigma value (needs reference)
            # self.peak_y_gussian = peak_height*(1/(hwhm_sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((self.x-peak_center)/hwhm_sigma)**2))) # Reference Gussian calculation
            peak_y_gussian = peak_height_a * np.exp(-((self.x - peak_center_b)**2)/(2 * (hwhm_sigma_c**2))) # Generic gaussian curve 
            self.all_gussian_y_arrays[i] = peak_y_gussian
            # peak_area_trapz = np.trapz(peak_y_gussian, self.x)
            # peak_area_anal = (math.sqrt(2) * math.pi) * (peak_height * hwhm_sigma)
            peak_area = (math.sqrt(2) * math.pi) * (peak_height_a * hwhm_sigma_c)
            # print(peak_area_trapz, peak_area_anal) 
            self.all_peak_areas[i] = peak_area
            peak_wavenumber = self.x[peak]
            self.all_peak_wavenumbers[i] = peak_wavenumber
            temp_array = np.array([peak_height_a, peak_center_b, hwhm_sigma_c])
            self.all_gussian_components[i] = temp_array

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
    
    # def peak_information(self):
        # spectrum_x_array = self.x
        # all_gaussian_y_arrays = self.all_gussian_y_arrays
        # peak_indecies = self.all_peaks
        # n = peak_indecies.size
        # self.all_peak_areas = np.empty(n)
        # self.all_peak_wavenumbers = np.empty(n)
        # for i in range(n):
            # peak = peak_indecies[i]
            # peak_y = all_gaussian_y_arrays[i]
            # peak_area = 
            # peak_area = np.trapz(peak_y, spectrum_x_array)
            # self.all_peak_areas[i] = peak_area
            # peak_wavenumber = self.x[peak]
            # elf.all_peak_wavenumbers[i] = peak_wavenumber

    def peak_comparison(self):
        # x_array for plotting curve
        self.x_curve_array = np.linspace(700, 1200, (self.x.size))
        # Importing internal spectrums 
        directory = "Internal Peak Database" # designating the project data folder 
        file_list = [os.path.join(directory, file) for file in os.listdir(directory)] # Array of all file paths
        n = len(file_list)
        self.all_curve_results = []
        import_spectrum_components = self.all_gussian_components
        for l in range(len(import_spectrum_components)):
            D, E, F = import_spectrum_components[l]
            self.all_curve_results.append(Curve_Results(Curve(D,E,F)))
        # looping through all internal spectrums
        for i in range(n):
            # Import the spectrum
            file_name_full = file_list[i]
            file_name = file_name_full.split('\\')[-1]
            internal_df = pd.read_csv(file_list[i], delim_whitespace=True, skiprows=1, header=None)
            internal_df.rename(columns={0: 'A', 1: 'B', 2: 'C'}, inplace=True)
            internal_spectrum_components = internal_df.to_numpy()
            comparison_array = np.empty(shape=[import_spectrum_components.shape[0], internal_spectrum_components.shape[0]])
            # Curve results object creation for all peaks in imported curve
            # looping through peak components
            for j in range(len(import_spectrum_components)):
                imp_curve = self.all_curve_results[j].curve
                # imp_curve_y = imp_curve.GetY(self.x_curve_array)
                for k in range(len(internal_spectrum_components)):
                    A, B, C = internal_spectrum_components[k]
                    int_curve = Curve(A, B, C)
                    # int_curve_y = int_curve.GetY(self.x_curve_array)
                    area = imp_curve.GetIntersectionArea(int_curve) # getting the intercetion area of the curves
                    comparison_array[j, k] = area
                    self.all_curve_results[j].internal_comp[file_name] = comparison_array[j]
        # Need to find the greatest intersection area for all values within one peak curve results
        for i in range(len(import_spectrum_components)):
            self.all_curve_results[i].FindHighestIntersectionArea()

    def spectrum_peak_output(self):
        gussian_peak_components = self.all_gussian_components
        output_df = pd.DataFrame(gussian_peak_components)
        output_df.columns = ["Peak Height (a)", "Peak Center (b)", "Peak Sigma (c)"]
        output_df.to_csv(f'Saved Output/Spectrum {self.file_index} Component Output.csv', sep='\t', index=False)

    def plot(self):
        # -- Whole spectrum figure
        fig = plt.figure(1)
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True)
        fig.suptitle('Full Spectrum Plot', fontsize=20, fontweight='bold')
        axs[0].set_title(self.file_name, fontsize=10)
        # -- Specific Plots
        axs[0].plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data') # Plotting of the spectrum
        for y_gussian_component in self.all_gussian_components:
            A, B, C = y_gussian_component
            axs[0].plot(self.x, GetY(self.x_curve_array, A, B, C), linewidth=1, linestyle='dashdot', color='black') # Plotting fitted gussian curves
        # for y2_gussian_array in sec_spectrum.all_gussian_y_arrays:
            # axs[0].plot(sec_spectrum.x, y2_gussian_array, linewidth=1, linestyle='dashdot', color='blue', label=f'{sec_spectrum.file_name} Peaks')
        # -- Plotting peak wavenumber\
        peak_number = self.all_peaks.size
        for peak in self.all_peaks:
            x, y = self.x[peak], self.y[peak]
            axs[0].text((x+3),(y-0.1), f'({peak_number})')
            peak_number -= 1
        # peak_number = sec_spectrum.all_peaks.size
        # for peak in sec_spectrum.all_peaks:
            # x, y = sec_spectrum.x[peak], sec_spectrum.y[peak]
            # axs[0].text((x+3),(y-0.1), f'({peak_number})')
            # peak_number -= 1
        if peak_deconvolution_type == 'first': # plotting all the first derv data 
            axs[1].plot(self.x, self.y_first_derivative, color='blue', linewidth=1.5, label='First Derivative') # plotting first dev on sep axs
            axs[1].plot(self.x, np.zeros_like(self.x), color='black', linewidth=1) # visualising zero line intersections
            axs[0].scatter(self.x[self.all_peaks], self.y[self.all_peaks], color='Black', marker=7, label='Peak Center')
            axs[1].scatter(self.x[self.first_dev_intersection_indices], self.y_first_derivative[self.first_dev_intersection_indices], color='black') 
            # # Visual of first derv peak locations
        if peak_deconvolution_type == 'second': # plotting all the second derv data
            axs[1].plot(self.x, self.y_second_derivative, color='green', linewidth=1.5, label="Second Derivative")# Plotting the second derv on sep axs
            #axs[1].scatter(self.x[self.spectrum_peaks], self.y_second_derivative[self.spectrum_peaks], color='black') # Visual of second derv peak locations
            axs[0].scatter(self.x[self.all_peaks], self.y[self.all_peaks], color='Black', marker=7, label='Peak Center') # Visual of second derv peak location
            axs[1].scatter(self.x[self.all_peaks], self.y_second_derivative[self.all_peaks], color='Black', marker=6) # Visual of second derv peak location
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
        fig2, axs2 = plt.subplots(3, 3)
        fig2.suptitle('Individual Peaks Plot', fontsize=20, fontweight='bold')
        n = 0
        # TODO Comment this mess
        peak_number = self.all_peaks.size
        for a in range(3):
            for b in range(3):
                if n >= self.all_peaks.size:
                    break
                else:
                    D1, E1, F1 = self.all_gussian_components[n]
                    peak = self.all_peaks[n]
                    axs2[a, b].plot(self.x, self.y, color='red', linewidth=1.5, label='XY Data')
                    # test = GetY(x, D, E, F)
                    axs2[a, b].plot(self.x, GetY(self.x_curve_array, D1, E1, F1), linewidth=1, color='blue', label='Fitted Peak')
                    axs2[a, b].scatter(self.x[self.all_peaks[n]], self.y[self.all_peaks[n]], color='black', marker=7, label='Peak Center')
                    axs2[a, b].set_xlim([np.max(self.x), np.min(self.x)])
                    axs2[a, b].set_ylim(0.00, 1.40)
                    axs2[a, b].title.set_text(f'Peak Number: {peak_number}')
                    axs2[a, b].text((self.x.max()-1), (1.10), f'Wavenumber (cm-1): {np.rint(self.all_peak_wavenumbers[n])}\n'
                    f'Peak Area: {np.round(self.all_peak_areas[n], 2)}\n' f'Reference Peak: {self.all_curve_results[n].internal_comp_highest_value[0]}')
                    peak_number -= 1
                    n += 1
        handles, labels = axs2[0, 0].get_legend_handles_labels()
        fig2.legend(handles, labels, loc='upper right')
        peak_plotting_t.stop()
        plt.show()
# Main Loop - Processing the spectrum
print("---------- START ----------")
# TODO: Add method to init all functions within the class (automate)
# Assigning spectrum object names
pri_spectrum = Spectrum('Import')
sec_spectrum = Spectrum('reference')
pri_spectrum.file_import_selection()
pri_spectrum.x_y_assignment()
pri_spectrum.x_axis_limiting()
pri_spectrum.normalise_data()
pri_spectrum.peak_detection()
pri_spectrum.peak_deconvolution()
pri_spectrum.peak_comparison()
# pri_spectrum.spectrum_peak_output()
# sec_spectrum.file_import_selection()
# sec_spectrum.x_y_assignment()
# sec_spectrum.x_axis_limiting()
# sec_spectrum.normalise_data()
# sec_spectrum.peak_detection()
# sec_spectrum.peak_deconvolution()
# sec_spectrum.peak_comparison()
pri_spectrum.plot()

total_runtime_t.stop()
print("-----------END-----------")