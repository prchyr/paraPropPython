import paraProp as prop
import matplotlib.pyplot as plt
import numpy as np
import signal
import sys



if len(sys.argv)<4:
    print("usage: python3 simpleExample.py <freq [GHz]> <depth [m]> <density fluctuations? [0,1]>")
    exit()
    
#set the simulation parameters
freq=float(sys.argv[1]) #frequency GHz
xMax=300 #m
zMax=200 #m
sourceDepth=float(sys.argv[2]) #m

polarization=0#vertical
site="SP"# south pole
method="II"# in ice
randomDensity=int(sys.argv[3])
nProfile="functional"#functional form for n(z)
if randomDensity==1:
    nProfile="data"#data-driven n(z) profile

#declare a simple parabolic solver
simple=prop.paraPropSimple(freq, xMax, zMax, sourceDepth, polarization, site, method, nProfile)


simple.doSolver()

#plot the magnitude of the field below the ice and slightly above
absu=(abs(simple.psiFull))
plt.style.use(['dark_background'])
plt.imshow(absu, extent=(0, xMax-1, -2*zMax, 2*zMax), aspect='auto', cmap='hot',  vmin=.0001, vmax=.02)
plt.title(str(int(freq*1000))+" MHz")
plt.xlabel("x [m]")
plt.ylabel("z [m]")
plt.ylim(-zMax, 100)

plt.show()


