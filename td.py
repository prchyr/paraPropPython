# paraPropPython time-dependent signal example use of paraPropPython.py notebook
# s. prohira, c. sbrocco

import paraPropPython as ppp
import numpy as np
import matplotlib.pyplot as plt
import util as util
import time

##### time-dependent example #####

### first, initialize an instance of paraProp by defining its dimensions and frequency of interest ###
iceDepth = 200. # m
iceLength = 100. # m
dx = 1 # m
dz = 0.05 # m

freq = 0.15

### it is useful to set the reference depth as the source depth when you only have one transmitter ###
sourceDepth = 30. # m
sim = ppp.paraProp(iceLength, iceDepth, dx, dz, refDepth=sourceDepth, airHeight=1)

def southpole(z):
    A=1.78
    B=-0.43
    C=-0.0132
    return A+B*np.exp(C*z)
sim.set_n(nFunc=southpole)

sim.set_dipole_source_profile(freq, sourceDepth)

### set a td signal ###
dt = 1
impulse = np.zeros(2**8, dtype='complex')
impulse[10] = 1+0j
sig = util.normToMax(util.butterBandpassFilter(impulse, 0.09, 0.25, 1/dt, 4))
sim.set_td_source_signal(sig, dt)

rxList = [ppp.receiver(100, 25)]
tic = time.perf_counter()
### run the solver ###
sim.do_solver(rxList)
toc = time.perf_counter()
print(f"Time: {toc-tic:0.4} seconds")

rx = rxList[0]
t = rx.get_time()
sig = rx.get_signal().real
f = rx.get_frequency()
spec = rx.get_spectrum()

wrapped = np.roll(sig, -np.argmax(sig)+45)
plt.plot(t, wrapped/max(wrapped))
plt.xlabel("Time (ns)")
plt.ylabel("Field (norm.)")
plt.show()
plt.plot(f, abs(spec)) 
plt.xlabel("Frequency (GHz)")
plt.ylabel("Mag. (abs)")
#plt.xlim(0,0.5)
plt.show()