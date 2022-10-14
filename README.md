# paraPropPython

this is a simple parabolic equation EM solver, it uses the parabolic equation approximation to simulate the propagation of EM waves in a medium with a changing index of refraction. 
Currently it is designed to simulate propagation of waves beneath the ice and in the air in the UHF frequency range on relatively short baselines 

This is a modified version of the original by Prohira and Cade Sbrocco. This version includes a number of updates including:
1. Two dimensional refractive index profiles (range varying)
2. Complex refractive index profiles (models wave attenuation)
3. Backwards reflected waves
4. Adds multi-core processing -> to allow more efficient computation

## Installation

no installation, just clone the repo

### Required packages
The code requires several python modules to run successfully:

These include:
* numpy
* scipy
* matplotlib
* shapely

The recommendation is to create a Anaconda or Miniconda environment. You can download and install Miniconda via: https://docs.conda.io/en/latest/miniconda.html



## using paraPropSimple

cd into the repo directory and try:

python3 simpleExample.py <frequency [GHz, keep it below 1]> <source depth [m]> <use density fluctiations? [0=no, 1=yes]>

and you should see a plot. 

## Time-domain / A-Scans
The simple version runs a simulation for a single frequency and a single source depth. 
Analysis of signal propagation through ice requires time-domain simulations. 
This is implemented by decomposing an input signal (i.e. a gaussian pulse) into a complex spectrum via Fourier transform.
The frequency dependent complex amplitudes are used to solve the EM field (for a given source dpeth) across the simulation geometry.
The field amplitude is sampled at some number of positions in the ice ('receivers') -> with each receiver having a received spectrum.
The receiver spectrum is recomposed via inverse Fourier transform into a receiver signal. The received signal is known as an 'A-Scan'

Am example of this is shown in example-td.ipynb

## B-Scans
Surface and borehole penetrating radar (SPR and BHR) is a tried and tested technique for non-invasive searches for interesting objects within the Earth's surface and the surfaces of other planets.
A common approach in SPR/BHR is Synthetic Aperture Radar (SAR) -> which utilizes a moving transceiver (a combined transmitting/TX and receiving antenna/RX).
In SPR, this is accomplished by moving the transceiver horizontally across the surface (i.e. in a straight line or a circle).
By combining a multitude A-Scans (measured time-domain waveform) for different transceiver positions, 
it's possible to create an image of the sub-surface, with targets being visible from the radar echo. These images are known as 'B-Scans'.
A BHR works by moving a transceiver vertically through a borehole -> creating a image of radar echos for different depths (this technique is ideal for identifying cracks and layers).

At present, a BHR can be simulated by running multiple PE simulations, where the source depth is moving vertically. For each source depth, E-field is across the whole spectrum, sampled at the receiver positions and recomposed into a receiver signal.
A BHR scan is created by combing the A-Scan of receivers at equal depth to the transmitter.

Simulating a SPR is not yet implemented.

TODO: Add an example of how to run a B-Scan

## paraPropBackwards
this is an update of the original code to include a range dependent refractive index (n = n(x,z)) and to allow simulation of backwards scattered EM waves off of objects by implementing a 2 way split step parabolic equation equation solver.

Simulation data, metadata and settings are haved to HDF files
Simulations are set up using a 'input' text file. See example.
To run a simulation:
python runSimulation.py <input-file.txt> <output-file.h5>

## paraPropPython Multi-Sim
This version allows for efficient computation by spreading PE solvers for different frequencies (and or source depths) to different CPU cores, allowing the different simulations to be solved in paralell.
