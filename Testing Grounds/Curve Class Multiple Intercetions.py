from os import name
import numpy as np
from numpy.core.numeric import cross
from scipy import stats
import math
import matplotlib.pyplot as plt

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
      # work that out here
      return -1

    elif (len(intersections) == 1):
      left = self if self.b <= other.b else other
      right = self if self.b > other.b else other
      left_area = left.GetDefiniteArea(intersections[0], 10000)
      right_area = right.GetDefiniteArea(-10000, intersections[0])
      return left_area + right_area

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


curve1 = Curve(1, 1100, 4) # Base Curve
curve2 = Curve(1, 1110, 4) # 1 intersection, to the right
curve3 = Curve(1, 1090, 4) # 1 intersection, to the left
curve4 = Curve(.5, 1100, 10) # 2 intersections
curve5 = Curve(.9, 1100, 4) # 0 intersections
curve6 = Curve(0.9994043855522292, 961.0, 14.438470604902145)
curve7 = Curve(0.8624232633279483, 956.7620112, 4.095737302880505)


toPlot1 = curve6
toPlot2 = curve7
print("toPlot1 area total = ", toPlot1.GetTotalArea())
print("toPlot2 area total = ", toPlot2.GetTotalArea())
print("Intersection Points = ", toPlot1.GetIntersectionPoints(toPlot2))
print("Intersection Area = ", toPlot1.GetIntersectionArea(toPlot2))

x = np.linspace(700, 1200, 1000)
plt.figure()
plt.plot(x, toPlot1.GetY(x), color="blue", label="y1")
plt.plot(x, toPlot2.GetY(x), color="red", label="y2")
plt.plot()
plt.vlines(toPlot1.GetIntersectionPoints(toPlot2), 0, 1, colors="Black")
plt.legend()
plt.show()