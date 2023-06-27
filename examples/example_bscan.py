import numpy as np
import matplotlib.pyplot as pl
import sys
import matplotlib
from scipy import signal

import h5py
import sys
import os
import time
import datetime
sys.path.append('../')
import util
import paraPropPython as ppp
from receiver import receiver as rx
from transmitter import tx_signal
from data import create_memmap, create_hdf_bscan

#One Way Bscan

if len(sys.argv) == 7:
    fname_profile = sys.argv[1]
    fname_output = sys.argv[2]
    freq_input = float(sys.argv[3])
    Band = float(sys.argv[4])
    dz_src = float(sys.argv[5])
    maxSourceDepth = float(sys.argv[6])
else:
    print('error! please enter td_simul.py <path/to/eps_file.txt> <path/to/output> <freq [GHz]> <Bandiwidth [GHz]> <dz [m]> [iceDepth m]')
    sys.exit()


def main(fname):
    #Set Refractive Index Profile
    profile_data = np.genfromtxt(fname)
    n_profile = profile_data[:,1]
    z_profile = profile_data[:,0]

    ### first, initialize an instance of paraProp by defining its dimensions and frequency of interest ###
    # Simulation Geometry Parameters
    #====================================================================
    iceDepth = 20.0
    iceLength = 50.  # m
    dx = 1  # m
    airHeight = 10
    dz = 0.02
    dz_src = 0.5
    #====================================================================


    #Siganl Parameters
    #====================================================================
    freq_centre = freq_input
    B = Band
    freq_base = freq_centre - B/2
    freq_top = freq_centre + B/2
    freq_HP = freq_centre - B
    freq_LP = freq_centre + B

    # Create Pulse

    dt = 0.2 #Sample Time Interval
    amplitude = 1.0 #amplitude
    t_max = iceLength * 2 * 1.8 / util.c_light #total time
    t_offset = 10. #time of maximum pulse peak
    sig_tx = tx_signal(amplitude=amplitude, frequency=freq_centre, bandwidth=B, t_centre=t_offset, tmax=t_max, dt=dt)
    signal_pulse = sig_tx.get_gausspulse()
    nSamples = sig_tx.nSamples
    #====================================================================

    #Set Receiver Array
    #====================================================================
    minSourceDepth = -1.5
    minRange = 2.0
    dRange = 2.0

    source_depths = np.arange(minSourceDepth, maxSourceDepth, dz_src) #Array containing source depths
    rx_x = np.arange(minRange, iceLength, dRange) #Array containing receiver ranges
    rx_z = np.arange(minSourceDepth, iceDepth, dz_src) #Array containing receiver depths
    nRX_x = len(rx_x)
    nRX_z = len(rx_z)
    nSources = len(source_depths)
    #====================================================================

    #Create Memmap and H5py file to save simulation results
    fname_npy = fname_output + '.npy'
    fname_h5 = fname_output + '.h5'

    #Memmap
    bscan_sig = create_memmap(fname_npy, dimensions=(nSources, nRX_x, nRX_z, nSamples))

    for i in range(nSources):
        sourceDepth = source_depths[i]
        sim = ppp.paraProp(iceLength=iceLength, iceDepth=iceDepth, dx=dx, dz=dz) #Set Simulation Geometry
        sim.set_n(nVec=n_profile, zVec=z_profile) #Set Refractive Index Profile
        sim.set_dipole_source_profile(freq_centre, sourceDepth)  # Set Source Profile
        sim.set_td_source_signal(signal_pulse, dt) #Set transmitted signal

        # Create HDF File -> Stores all simulation data and metadata (must do this on the first iteration and only on the first iteration)
        if i == 0:
            output_hdf = create_hdf_bscan(fname_h5, sim=sim, tx_signal=sig_tx, tx_depths=source_depths, rx_ranges=rx_x,
                                      rx_depths=rx_z)

        #Create List of Receivers
        rxList0 = []
        for j in range(nRX_x):
            for k in range(nRX_z):
                rx_jk = rx(x=rx_x[j], z=rx_z[k])
                rxList0.append(rx_jk)

        print('Start solution')
        tstart = time.time() #Records Simulation Time

        #Run Simulation
        sim.do_solver(rxList=rxList0, freqMin=freq_HP, freqMax=freq_LP)
        #Solution Done
        tend = time.time()

        ii = 0

        #Iterates over each receiver in list -> saves receiver signal to the memmap [transmitter_id, rx_range_id, rx_depth_id]
        for j in range(nRX_x):
            for k in range(nRX_z):
                rx_jk = rxList0[ii]
                sig_rx = rx_jk.get_signal()
                nSig = len(sig_rx)
                if nSamples != nSig:
                    dN = nSamples - nSig
                    if dN > 0:
                        bscan_sig[i, j, k, :nSig] = sig_rx
                    else:
                        bscan_sig[i, j, k, :] = sig_rx[:nSamples]
                else:
                    bscan_sig[i, j, k, :] = sig_rx
                ii += 1
        duration = tend - tstart
        print('Solution complete, duration = ', datetime.timedelta(seconds=duration))

        remainder = (nSources - i) * duration
        print('Remaining scan time: ', datetime.timedelta(seconds=remainder))
        print('')
    #Data is added to HDF file
    output_hdf.create_dataset('bscan_sig', data=bscan_sig)
    output_hdf.close()
    print('Bscan complete')

    return -1


if __name__ == "__main__":
    main(fname_profile)