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

path2jobs = fname_out[:-3]
if os.path.isdir(path2jobs) == False:
    os.system('mkdir ' + path2jobs)
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

fname_matrix = fname_out[:-3] + '_matrix.h5'
output_hdf = h5py.File(path2jobs + '/' + fname_matrix,'w')
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


fname_joblist = path2jobs + '/' + fname_out[:-3] + '_A=' + str(A_parameter) + '_joblist.txt'
job_list = open(fname_joblist,'w')

#Run Reference File
fname_ref = fname_out[:-3] + '-ref.h5'
fname_ref_with_path = path2jobs + '/' + fname_ref
if os.path.isfile(fname_ref_with_path) == False:
    create_ref_profile = 'python runSim_for_nprofile_ref.py ' + fname_config + ' ' + fname_profile + ' ' + fname_ref_with_path
    os.system(create_ref_profile)

for i in range(nRandom):
    jobline = 'python runSim_for_nprofile.py ' + fname_ref + ' ' + fname_matrix + ' ' + str(i) + '\n'
    job_list.write(jobline)
job_list.close()