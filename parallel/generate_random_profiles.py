import numpy as np
import matplotlib.pyplot as pl
import os
from scipy.interpolate import interp1d
import datetime
import time
import sys
import h5py


sys.path.append('../')
import paraPropPython as ppp
from transmitter import tx_signal
from receiver import receiver as rx
from data import create_memmap, create_hdf_bscan
from permittivity import southpole
from permittivity import rho2n
import util
import configparser

B = 1.0
C = 0.01
D = 0.5
E = 1.0
low_cut = 0.5

def makeRandomDensityVector(z, a=0.6, b=B, c=C, d=D, e=E, low_cut=low_cut):
    """make a vector of random density fluctuations. This is currently used with the Taylor Dome n(z) profile.
    the density fluctuations are on the order of a few percent of the density."""
    dz = abs(z[1] - z[0])
    ranVec = util.lowpassFilter(dz, low_cut, (a / (b + (z * c))) * (e*np.random.random_sample(len(z)) - d))
    return ranVec

if len(sys.argv) != 6:
    print('error, enter python ' + sys.argv[0] + ' <nprofile.txt> <config.txt> <fname_out.h5> <A> <nRandom>')
    sys.exit(-1)
fname_profile = sys.argv[1]
fname_config = sys.argv[2]
fname_out = sys.argv[3]
A_parameter = float(sys.argv[4])
nRandom = int(sys.argv[5])

n_data = np.genfromtxt(fname_profile)
z_profile = n_data[:,0]
n_profile_0 = n_data[:,1]
nDepths = len(z_profile)

config = configparser.ConfigParser()
config.read(fname_config)
geometry = config['GEOMETRY']
iceDepth = float(geometry['iceDepth'])
ii_cut = util.findNearest(z_profile, iceDepth)
nDepths_out = ii_cut
n_random_array = np.ones((nRandom, nDepths_out),dtype='float')
z_profile_out = np.linspace(0, iceDepth, nDepths_out)

n_random_array[0] = z_profile_out
for i in range(1,nRandom):
    ranVec = makeRandomDensityVector(z_profile, a=A_parameter)
    n_profile_i = n_profile_0 + ranVec
    n_profile_out = n_profile_i[:ii_cut]
    n_random_array[i,:] = n_profile_out

output_hdf = h5py.File(fname_out,'w')
output_hdf.attrs['iceDepth'] = iceDepth
output_hdf.attrs['A'] = A_parameter
output_hdf.attrs['B'] = B
output_hdf.attrs['C'] = C
output_hdf.attrs['D'] = D
output_hdf.attrs['E'] = E
output_hdf.attrs['low_cut'] = low_cut

output_hdf.attrs['nRandom'] = nRandom
output_hdf.attrs['n profile filename'] = fname_profile
output_hdf.create_dataset('n_profile_random',data=n_random_array)
output_hdf.create_dataset('n_profile_input', data=n_profile_0)
output_hdf.create_dataset('z_profile_input', data=z_profile)
output_hdf.create_dataset('z_profile_output', data=z_profile_out)

S_arr = np.zeros(nRandom)
output_hdf.create_dataset('S_arr',data=S_arr)
output_hdf.close()