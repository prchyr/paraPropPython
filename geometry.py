import numpy as np
from math import pi
from permittivity import eps2m
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#TODO Create a general geometry class -> create closed/polygon geometry from a list of coordinates
#TODO Create a irregular surface class -> create surface geometry from a list of coordiantes

def det(v1, v2):
    V1 = [v1.dx, v1.dz]
    V2 = [v2.dx, v2.dz]
    return np.cross(V1, V2)

def area(p1, p2, p3):
    A = abs((p1.x * (p2.z - p3.z) + p2.x * (p3.z - p1.z) + p3.x * (p1.z - p2.z)) / 2.0)
    return A

class point:
    def __init__(self, x, z):
        self.x = x
        self.z = z

class vector:
    def __init__(self, p0, p1):
        self.x0 = p0.x
        self.z0 = p0.z
        self.x1 = p1.x
        self.z1 = p1.z
        self.dx = self.x1 - self.x0
        self.dz = self.z1 - self.z0
        self.mag = np.sqrt(self.dx ** 2 + self.dz ** 2)
        self.arg = np.arctan2(self.dz, self.dx)

class circle:
    def __init__(self, x, z, r, eps_r = 1+0j):
        self.x = x  # X_position central
        self.z = z  # Z_position central
        self.r = r  # Radius
        self.eps_r = eps_r
        self.n = eps2m(eps_r)
        self.A = 4 * pi * self.r ** 2
        self.V = (4. / 3.) * pi * self.r ** 3

    def isInside(self, x1, z1):
        dx = x1 - self.x
        dz = z1 - self.z
        dr = np.sqrt(dx ** 2 + dz ** 2)
        if dr < self.r:
            return True
        else:
            return False

class triangle:
    def __init__(self, p1, p2, p3, eps_r = 1+0j):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.area = area(p1, p2, p3)

        self.eps_r = eps_r
        self.n = eps2m(self.eps_r)

    def isInside(self, x, z):
        p = point(x,z)
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3

        # Calculate area of triangle ABC
        A = area(p1, p2, p3)

        # Calculate area of triangle PBC
        A1 = area(p, p2, p3)

        # Calculate area of triangle PAC
        A2 = area(p1, p, p3)

        # Calculate area of triangle PAB
        A3 = area(p1, p2, p)

        # Check if sum of A1, A2 and A3
        # is same as A
        #print(A, A1, A2, A3)
        if (A == A1 + A2 + A3):
            return True
        else:
            return False

class polygon:
    def __init__(self, list_points, eps_r=1+0j):
        self.npoints = len(list_points)
        self.eps_r = eps_r
        self.n = eps2m(eps_r)
        self.shape = Polygon(list_points)

    def isInside(self, x, z):
        pt = Point(x,z)
        return self.shape.contains(pt)

