import sys
import numpy as np
import time
import datetime

sys.path.append('../')
import paraPropPython as ppp
from receiver import receiver as rx
from transmitter import tx_signal
from data import create_sim, create_rx_ranges, create_hdf_bscan
from data import create_transmitter_array, create_tx_signal, create_memmap, create_rxList

#I want -> runSum_for_nprofile.py <config.txt> <nprofile.txt> <output.txt>

if len(sys.argv) != 4:
    print('error! you must enter argument: \npython runSim_for_nprofile.py <fname_config.txt> '
          '<fname_nprofile.txt>')

fname_config = sys.argv[1] #cofniguration file
fname_nprofile = sys.argv[2] # n-profile.txt
fname_out = sys.argv[3] # output-file.h5

tx_signal = create_tx_signal(fname_config)
tx_signal.get_gausspulse()

rx_ranges = create_rx_ranges(fname_config)
tx_depths = create_transmitter_array(fname_config)
nDepths = len(tx_depths)
rx_depths = tx_depths
nRX_x = len(rx_ranges)
nRX_z = len(rx_depths)

n_data = np.genfromtxt(fname_nprofile)
z_profile = n_data[:,0]
n_profile = n_data[:,1]

print(tx_signal.freqMin, tx_signal.freqMax)

bscan_npy = np.zeros((nDepths, nRX_x, nRX_z, tx_signal.nSamples),dtype='complex')
for i in range(nDepths):

    tstart = time.time()
    #=========================================================================================================
    sourceDepth = tx_depths[i]
    print('start solution for ', fname_nprofile, ' z_tx =', sourceDepth, ' m below surface', i+1, ' out of', nDepths)

    sim = create_sim(fname_config)
    sim.set_n(nVec=n_profile, zVec=z_profile) #Set Refractive Index Profile
    sim.set_dipole_source_profile(tx_signal.frequency, sourceDepth)  # Set Source Profile
    sim.set_td_source_signal(tx_signal.pulse, tx_signal.dt) #Set transmitted signal
    rxList = create_rxList(rx_ranges, rx_depths)

    sim.do_solver(rxList, freqMin=tx_signal.freqMin, freqMax=tx_signal.freqMax)
    if i == 0:
        output_hdf = create_hdf_bscan(fname=fname_out, sim=sim, tx_signal=tx_signal,
                                      tx_depths=tx_depths, rx_ranges=rx_ranges, rx_depths=rx_depths)

    ii = 0
    for j in range(nRX_x):
        for k in range(nRX_z):
            rx_jk = rxList[ii]
            bscan_npy[i,j,k,:] = rx_jk.get_signal()
            ii += 1
    #==========================================================================================================
    tend = time.time()
    duration_s = (tend - tstart)
    duration = datetime.timedelta(seconds=duration_s)
    remainder = duration_s * (nDepths - (i+1))
    completed = round(float(i+1)/float(nDepths) * 100,2)
    print(completed, ' % completed')
    print('remaining steps', nDepths - (i+1),'\nremaining time:', remainder,'\n')
    now = datetime.datetime.now()
    tstamp_now = now.timestamp()
    end_time = datetime.datetime.fromtimestamp(tstamp_now + duration_s)
    print('completion at:', end_time)
    print('')
output_hdf.create_dataset('bscan_sig', data=bscan_npy)
output_hdf.close()

