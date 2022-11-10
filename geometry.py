import numpy as np
from math import pi
from permittivity import eps2m
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#TODO Create a irregular surface class -> create surface geometry from a list of coordiantes

'''
Geometry Module:
This includes functions for creating 2D refracitve index profiles
'''

def det(v1, v2):
    '''
    Calculates the cross product between two 2D-vectors
    #TODO: Consider renmaing this -> name is misleading
    '''
    V1 = [v1.dx, v1.dz]
    V2 = [v2.dx, v2.dz]
    return np.cross(V1, V2)

def area(p1, p2, p3):
    '''
    Calculates the area of a triangle with 3 points all contained in a 2D surface
    p1, p2, p3 : Points with 2 coordinates (x,z)
    '''
    A = abs((p1.x * (p2.z - p3.z) + p2.x * (p3.z - p1.z) + p3.x * (p1.z - p2.z)) / 2.0)
    return A

class point:
    '''
    Defines a point on a 2D surface
    '''
    def __init__(self, x, z):
        self.x = x
        self.z = z

class vector:
    '''
    Defines a 2D vector from 2 points
    '''
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
    '''
    Defines a circle -> used to define a circular object with a constant permittivity
    i.e. can represent the cross section of a pipe or a sphere in a 2D permittivity profile

    Inputs:
    x, z : defines centre of circle
    r : radius of circle
    eps_r : permittivity of circular object
    '''
    def __init__(self, x, z, r, eps_r = 1+0j):
        self.x = x  # X_position central
        self.z = z  # Z_position central
        self.r = r  # Radius
        self.eps_r = eps_r
        self.n = eps2m(eps_r)
        self.A = 4 * pi * self.r ** 2
        self.V = (4. / 3.) * pi * self.r ** 3

    def isInside(self, x1, z1):
        '''
        Checks in a point defined by x1 and z1 is inside of the circle
        '''
        dx = x1 - self.x
        dz = z1 - self.z
        dr = np.sqrt(dx ** 2 + dz ** 2)
        if dr < self.r:
            return True
        else:
            return False

class triangle:
    '''
        Used to define a triangular object with a constant permittivity
        i.e. can represent a canyon, crevass or water-filled crack

        Inputs:
        p1,p2,p3 -> points defining the triangle
        eps_r : permittivity of circular object
        '''
    def __init__(self, p1, p2, p3, eps_r = 1+0j):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.area = area(p1, p2, p3)

        self.eps_r = eps_r
        self.n = eps2m(self.eps_r)

    def isInside(self, x, z):
        '''
        Checks in a point defined by x and z are inside of the circle
        '''
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
    '''
    Defines an arbitary 2D shape with constant permittivity
    -> can represent the cross-section of a volume of constant internal permittivity

    inputs:
    list_points : a python list containing points defining the premitter of the object (should be a closed path)
    eps_r : internal refractive index

    '''
    def __init__(self, list_points, eps_r=1+0j):
        self.npoints = len(list_points)
        self.eps_r = eps_r
        self.n = eps2m(eps_r)
        self.shape = Polygon(list_points)

    def isInside(self, x, z):
        pt = Point(x,z)
        return self.shape.contains(pt)

