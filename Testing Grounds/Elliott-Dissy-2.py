from os import name
import numpy as np
from numpy.core.numeric import cross
from scipy import stats
import math
import matplotlib.pyplot as plt
import pandas as pd

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

# ---------- Importing curve data from database:
Meteorite_Curves_df = pd.read_csv('Saved Output\Spectrum 12 Component Output.csv', delim_whitespace=True, skiprows=1, header=None)
Meteorite_Curves_df.rename(columns={0: 'a', 1: 'b', 2: 'c'}, inplace=True)
Import_Curve_Data = Meteorite_Curves_df.to_numpy()
Augite_Curves_df = pd.read_csv('Saved Output\Spectrum 13 Component Output.csv', delim_whitespace=True, skiprows=1, header=None)
Augite_Curves_df.rename(columns={0: 'a', 1: 'b', 2: 'c'}, inplace=True)
Reference_Curve_Data = Augite_Curves_df.to_numpy()

# ---------- Comparing one peak against all peaks:
n = 0
x = np.linspace(700, 1200, 1000)
for a in Import_Curve_Data: # Create curve info for import peak:
  ahhhhhhh = Reference_Curve_Data.shape[0]
  fig, axs = plt.subplots((ahhhhhhh), sharex=True)
  fig.suptitle('Individual Peaks Plot', fontsize=20, fontweight='bold')
  if n == ahhhhhhh:
    n = 0
  Imp_Curve = Curve(a[0], a[1], a[2])
  toPlot1 = Imp_Curve
  print("toPlot1 area total = ", toPlot1.GetTotalArea())
  imp_curve_y = toPlot1.GetY(x)
  one_half_sigma1 = toPlot1.c * 2.5
  for b in Reference_Curve_Data:
    print(f'Peak {n}')
    Ref_Curve = Curve(b[0], b[1], b[2])
    toPlot2 = Ref_Curve
    ref_curve_y = toPlot2.GetY(x)
    one_half_sigma2 = toPlot2.c * 2.5
    print("toPlot2 area total = ", toPlot2.GetTotalArea())
    print("Intersection Points = ", toPlot1.GetIntersectionPoints(toPlot2))
    print("Intersection Area = ", toPlot1.GetIntersectionArea(toPlot2))
    axs[n].plot(x, ref_curve_y, color='blue', label='Reference Spectrum')
    axs[n].plot(x, imp_curve_y, color='red', label='Import Spectrum')
    axs[n].vlines(toPlot1.GetIntersectionPoints(toPlot2), 0, 1, colors="Black", label='Intersection Points')
    axs[n].hlines(toPlot1.a, toPlot1.b - one_half_sigma1, toPlot1.b + one_half_sigma1)
    axs[n].hlines(toPlot2.a, toPlot2.b - one_half_sigma2, toPlot2.b + one_half_sigma2)
    n += 1
  handles, labels = axs[0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right')
  plt.show()