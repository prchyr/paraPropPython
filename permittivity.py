from numpy import sqrt
from math import pi
import numpy as np
import h5py
import util

#Fundamental Constants
epsilon0 = 8.85e-12 #Permittivity of Free Space
mu0 = 4*pi*1e-7 #Permeability of Free Space
c0 = 1./sqrt(epsilon0*mu0) #Speed of Light in Vacuum

#Functions
def eps2m(eps_c): #Return complex refractive index (m) from complex relative permittivity (eps_c)
    eps_r = eps_c.real
    eps_i = eps_c.imag
    n_sq = (1/2.) * (sqrt(eps_r**2 + eps_i**2) + eps_r )
    k_sq = (1/2.) * (sqrt(eps_r**2 + eps_i**2) - eps_r )
    n = sqrt(n_sq)
    k = sqrt(k_sq)
    m = n + 1j*k
    return m

def m2eps(m):
    eps_c = m**2
    return eps_c

def cond2eps_im(sigma, f):
    w = 2*pi*f
    return sigma/(w*epsilon0)

def sigma_from_p(p, sigma): #p: porosity, sigma: conductivity for p = 0
    v = 1 - p #volume filling factor
    return sigma*v*(0.68 + 0.32*v)**2

def alpha(eps_c, f):
    m_c = eps2m(eps_c)
    k = m_c.imag
    w = 2*pi*f
    return 4*w*k/c0