from os import name
import numpy as np
from numpy.core.numeric import cross
from scipy import special
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
    left = self if self.b <= other.b else other
    right = self if self.b > other.b else other
    # find the right intersection point
    intersection = left.GetIntersectionPoints(right)[1] # INTERRIM EXAMPLE!!!!!!!

    left_area = left.GetDefiniteArea(intersection, 10000)
    right_area = right.GetDefiniteArea(-10000, intersection)
    return left_area + right_area

  def GetIntersectionPoints(self, other):
    x = 1 / (2 * self.c ** 2) - 1 / (2 * other.c ** 2)
    y = other.b / (other.c**2) - self.b/(self.c**2)
    z = self.b**2 /(2*self.c**2) - other.b**2 / (2*other.c**2) - np.log((other.c*self.a)/(self.c*other.a))
    return np.roots([x,y,z])
  
  def GetY(self, x): # works for one point
    return self.a * np.exp(-((x - self.b) ** 2)/(2 * self.c ** 2))


curve1 = Curve(0.6969, 1075, 4.246)
curve2 = Curve(0.6809, 1095, 3.276)

print("Peak 1 area total = ", curve1.GetTotalArea())
print("Peak 2 area total = ", curve2.GetTotalArea())
print("Peak Intersection Area = ", curve1.GetIntersectionArea(curve2))
print("Peak Area Percentage Overlap = ", (curve1.GetIntersectionArea(curve2)/(curve1.GetTotalArea() + curve2.GetTotalArea())) * 100)


x = np.linspace(1000, 1200, 10000)
plt.figure()
plt.plot(x, curve1.GetY(x), color="blue", label="Peak One")
plt.plot(x, curve2.GetY(x), color="red", label="Peak Two")
plt.legend()
plt.show()



# cum = Curve(1, 2, 3)

# def solve(a1, a2, b1, b2, c1, c2):
#   x = 1/(2*c1**2) - 1/(2*c2**2)
#   y = b2/(c2**2) - b1/(c1**2)
#   z = b1**2 /(2*c1**2) - b2**2 / (2*c2**2) - np.log((c2*a1)/(c1*a2))
#   return np.roots([x,y,z])


# a1 = 0.6969
# b1 = 1075
# c1 = 4.246

# a2 = 0.6809
# b2 = 1095
# c2 = 3.276
# a=0
# b=0
# c=0

# secretanswer = -math.sqrt(math.pi) * math.sqrt(2) * (a2 * c2 - a1 * c1)

# AreaSharedByTwoGaussians = math.sqrt(math.pi) * math.sqrt(2) * (a2 * c2 - a1 * c1)
# AreaOfASingleGaussian = math.sqrt(math.pi) * math.sqrt(2) * (a1 * c1)

# print(AreaOfASingleGaussian)


# crossover = solve(a1, a2, b1, b2, c1, c2)
# print("crossover 1 -> 2 = ", solve(a1, a2, b1, b2, c1, c2))
# print("crossover 2 -> 1 = ", solve(a2, a1, b2, b1, c2, c1))

# x = np.linspace(1000, 1200, 10000)
# y1 = a1 * np.exp(-((x - b1) ** 2)/(2 * c1 ** 2))
# y2 = a2 * np.exp(-((x - b2) ** 2)/(2 * c2 ** 2))

# def integral(a, b, c, xmin, xmax):
#     root2 = math.sqrt(2)
#     p1 = math.erf((xmax - b) / (root2 * c))
#     p2 = math.erf((xmin - b) / (root2 * c))
#     return (math.sqrt(math.pi) / root2) * a * c * (p1 - p2)

# print("crossovers = ", crossover)
# intersection = crossover[1]

# arealeft = integral(a1, b1, c1, intersection, 1200)
# arearight = integral(a2, b2, c2, 1000, intersection)
# area_g1 = integral(a1, b1, c1, -10000, 10000)
# area_g2 = integral(a2, b2, c2, -10000, 10000)

# print("area of G1 (new) = ", area_g1)
# print("area of G2 (new) = ", area_g2)
# print("area intersection = ", arealeft + arearight)
# print("percentage = ", 100 * (arealeft + arearight) / (area_g1 + area_g2) )

# print("trapz y1 total = ", np.trapz(y1, x))
# print("trapz y2 total = ", np.trapz(y2, x))

# plt.figure()
# plt.plot(x, y1, color="blue", label="y1")
# plt.plot(x, y2, color="red", label="y2")
# plt.vlines(intersection, 0, 0.5)
# plt.legend()
# plt.show()


