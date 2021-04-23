import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import signal
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import itertools
# linearly spaced x-axis of 10 values between 1 and 10
x = np.linspace(1,100,1000)
x1 = x + 5

amp1 = 1000
amp2 = 2000
sigma1 = 5
cen1 = 50
wight_factor = 1
y_array_gauss = amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
y_array_gauss1 = amp2*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
# creating some noise to add the the y-axis data
y_noise_gauss = (np.exp((np.random.ranf(50))))/5
y_array_gauss1d = np.gradient(y_array_gauss)
y_array_gauss2d = np.gradient(y_array_gauss1d)
y_array_gauss1d = y_array_gauss1d * wight_factor
y_array_gauss2d = y_array_gauss2d * wight_factor
# peak testing:
def find_peak_maxima(x,y):
    peaks, _ = signal.find_peaks(y)
    x_peak = x[peaks]
    y_peak = y[peaks]
    return(peaks, x_peak, y_peak)

def find_peak_minima(x,y):
    y2 = y * -1
    peaks, _ = signal.find_peaks(y2)
    x_peak = x[peaks]
    y_peak = y[peaks]
    return(peaks, x_peak, y_peak)

def peak_height (y, peaks):
    ymax = y[peaks]
    ymin = y.min()
    return(ymax, ymin)


peaks_d0, x_d0, y_d0 = find_peak_maxima(x, y_array_gauss)
ymax_d0, ymin_d0 = peak_height(y_array_gauss, peaks_d0)
peaks_d1, x_d1, y_d1 = find_peak_maxima(x, y_array_gauss1d)
peaks_d2, x_d2, y_d2 = find_peak_maxima(x, y_array_gauss2d)
peaks_neg_d1, xneg_d1, yneg_d1 = find_peak_minima(x, y_array_gauss1d)


plt.plot(x, y_array_gauss, "r")
# plt.plot(x1, y_array_gauss1, "b")
# plt.fill_between(x, y_array_gauss, 0, color='lightcoral')
#plt.fill_between(x1, y_array_gauss, 0, color='lightblue')
# plt.fill_betweenx(y_array_gauss1, x, x1, where=(x!=x1), color='white', interpolate=True)
plt.plot(x, y_array_gauss1d, 'b')
plt.plot(x, y_array_gauss2d, 'g')
plt.scatter(x_d0, y_d0, color='black', s=101)
plt.scatter(x[peaks_d0], ymax_d0/2, color='black', s=101)
plt.scatter(x[peaks_d0], ymin_d0, color = 'black', s=101)
plt.scatter(x_d1, y_d1, color='black')
plt.scatter(xneg_d1, yneg_d1, color='black')
plt.vlines(x_d0, ymax_d0, ymin_d0, colors='black')
for peak in peaks_d2:
    plt.vlines(x[peak], y_array_gauss2d[peak], y_array_gauss[peak])
for peak in peaks_d1:
    plt.vlines(x[peak], y_array_gauss1d[peak], y_array_gauss[peak])
for peak in peaks_neg_d1:
    plt.vlines(x[peak], y_array_gauss1d[peak], y_array_gauss[peak])
plt.scatter(x_d2, y_d2, color='black')
plt.xlim(0, 100)
plt.show()