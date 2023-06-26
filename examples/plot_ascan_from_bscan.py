import numpy as np
import sys
sys.path.append('../')

import paraPropPython as ppp
from receiver import receiver
from data import bscan

import scipy.signal as signal
import util
import h5py
import matplotlib.pyplot as pl
import peakutils as pku
import sys

if len(sys.argv) == 5:
    fname = sys.argv[1]
    zTx = float(sys.argv[2])
    xRx = float(sys.argv[3])
    zRx = float(sys.argv[4])
    plot_option = 'real'
elif len(sys.argv) == 6:
    fname = sys.argv[1]
    zTx = float(sys.argv[2])
    xRx = float(sys.argv[3])
    zRx = float(sys.argv[4])
    plot_option = sys.argv[5]
else:
    print('error: must enter 4 or 5 arguments \n python plot_td_simul.py fname zTx xRx zRx <plot_option?>')
    sys.exit()

if plot_option != 'real' and plot_option != 'abs' and plot_option != 'dB':
    print('warning: plot options must be real, abs or dB \n default option: real')

fname_h5 = fname
sim_bscan = bscan()
sim_bscan.load_sim(fname_h5)
ascan = sim_bscan.get_ascan(zTx=zTx, xRx=xRx, zRx=zRx)
tspace = sim_bscan.tspace

pl.figure(figsize=(10,6),dpi=150)
pl.title('R = ' + str(round(xRx,2)) + ' m, z_tx = ' + str(round(zRx,2)) + ' m, z_rx = ' + str(round(zRx,2)))
if plot_option == 'real':
    pl.plot(tspace, ascan.real)
    pl.ylabel('Amplitude [V]')

elif plot_option == 'abs':
    pl.plot(tspace, abs(ascan)**2)
    pl.ylabel('Power [u]')
elif plot_option == 'dB':
    pl.plot(tspace, 20*np.log10(abs(ascan)))
    pl.ylabel('Power [dBu]')
else:
    pl.plot(tspace, ascan.real)
    pl.ylabel('Amplitude [V]')

n_profile = sim_bscan.n_profile
z_profile = sim_bscan.z_profile
tCentral = sim_bscan.tx_signal.t_centre
ii = util.findNearest(z_profile, zTx)
R = np.sqrt((zRx-zTx)**2 + xRx**2)
n_at_depth = n_profile[ii]
pl.axvline(tCentral + R*n_at_depth/util.c_light,c='k')
pl.xlabel('Time [ns]')
pl.grid()
pl.show()