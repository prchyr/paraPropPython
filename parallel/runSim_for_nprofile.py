import sys
import numpy as np
from matplotlib import pyplot as pl
import time
import datetime
import h5py

sys.path.append('../')
import paraPropPython as ppp
from receiver import receiver as rx
from transmitter import tx_signal
from data import create_sim, create_rx_ranges, create_hdf_bscan
from data import create_transmitter_array, create_tx_signal, create_memmap, create_rxList
from data import bscan
#I want -> runSum_for_nprofile.py <config.txt> <nprofile.txt> <output.txt>

if len(sys.argv) != 4:
    print('error! you must enter argument: \npython runSim_for_nprofile.py <reference_bscan.h5> '
          '<fname_nprofile_matrix.h5 i')

fname_reference = sys.argv[1] # comparison file
fname_n_matrix = sys.argv[2]
ii_select = int(sys.argv[3])

n_matrix_hdf = h5py.File(fname_n_matrix,'r')
n_profile_matrix = np.array(n_matrix_hdf.get('n_profile_random'))
n_profile_rand = n_profile_matrix[ii_select]
z_profile_rand = np.array(n_matrix_hdf.get('z_profile_output'))

n_matrix_hdf.close()

if ii_select >= len(n_profile_matrix):
    print('error! ii_select must be greater than zero and less than the number of randomized profiles in ', fname_n_matrix)
    print(0, ' < ii_select < ', len(n_profile_matrix))
    sys.exit(-1)
#fname_out = sys.argv[4]

bscan_ref = bscan()
bscan_ref.load_sim(fname_reference)
rx_ranges = bscan_ref.rx_ranges
rx_depths = bscan_ref.rx_depths
tx_signal = bscan_ref.tx_signal
tx_depths = bscan_ref.tx_depths
nRX_x = len(rx_ranges)
nRX_z = len(rx_depths)
nDepths = len(tx_depths)

z_profile_ref = bscan_ref.z_profile
n_profile_ref = bscan_ref.n_profile



bscan_npy = np.zeros((nDepths, nRX_x, nRX_z, tx_signal.nSamples),dtype='complex')
for i in range(nDepths):
    tstart = time.time()
    # =========================================================================================================
    sourceDepth = tx_depths[i]
    print('start solution for z_tx =', sourceDepth, ' m below surface', i+1, ' out of', nDepths)

    sim = ppp.paraProp(iceDepth=bscan_ref.iceDepth, iceLength=bscan_ref.iceLength, dx=bscan_ref.dx, dz=bscan_ref.dz, airHeight=bscan_ref.airHeight)
    sim.set_n(nVec=n_profile_rand, zVec=z_profile_rand) #Set Refractive Index Profile
    sim.set_dipole_source_profile(tx_signal.frequency, sourceDepth)  # Set Source Profile
    sim.set_td_source_signal(tx_signal.pulse, tx_signal.dt) #Set transmitted signal
    rxList = create_rxList(rx_ranges, rx_depths)

    sim.do_solver(rxList, freqMin=tx_signal.freqMin, freqMax=tx_signal.freqMax)

    ii = 0
    for j in range(nRX_x):
        for k in range(nRX_z):
            rx_jk = rxList[ii]
            bscan_npy[i,j,k,:] = rx_jk.get_signal()
            ii += 1
    # ==========================================================================================================
    tend = time.time()
    duration_s = tend - tstart
    duration = datetime.timedelta(seconds=duration_s)
    remainder = datetime.timedelta(seconds=duration_s * (nDepths - (i + 1)))
    completed = round(float(i + 1) / float(nDepths) * 100, 2)
    print(completed, ' % completed, solution time: ', duration)
    print('remaining steps', nDepths - (i + 1), '\nremaining time:', remainder, '\n')
    now = datetime.datetime.now()
    tstamp_now = now.timestamp()
    end_time = datetime.datetime.fromtimestamp(tstamp_now + duration_s)
    print('completion at:', end_time)
    print('')



def signal_correlation(sig_rx_ref, sig_rx_sim):
    corr_arr = sig_rx_ref * sig_rx_sim
    corr = sum(corr_arr)
    return corr

Corr = 0
for i in range(nDepths):
    for j in range(nRX_x):
        for k in range(nRX_z):
            Corr += signal_correlation(abs(bscan_ref.bscan_sig[i,j,k]), abs(bscan_npy[i,j,k]))
S_corr = 1/Corr
print(Corr, S_corr)
n_matrix_hdf = h5py.File(fname_n_matrix,'r+')
S_arr = n_matrix_hdf['S_arr']
S_arr[i] = S_corr
n_matrix_hdf.close()