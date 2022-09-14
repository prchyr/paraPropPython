import util
import math as m
import numpy as np
from inspect import signature
from scipy.interpolate import interp1d
from scipy import signal
# A. Kyriacou

class tx_signal:
    def __init__(self, frequency, bandwidth, t_centre, dt, tmax, amplitude = 1, freqMin=None, freqMax=None):
        self.amplitude = amplitude
        self.frequency = frequency
        self.bandwidth = bandwidth
        self.t_centre = t_centre
        self.dt = dt
        self.fsample = 1/dt
        self.freq_nyq = 1/(2*dt)
        self.tmax = tmax
        self.nSamples = int(tmax/dt)
        self.tspace = np.linspace(0, tmax, self.nSamples)
        self.freq_space = np.fft.fftfreq(self.nSamples, self.dt)
        if freqMax == None:
            self.freqMax = self.frequency - self.bandwidth/2
        else:
            self.freqMax = freqMax
        if freqMin == None:
            self.freqMin = self.frequency - self.bandwidth/2
        else:
            self.freqMin = freqMin

    def get_gausspulse(self, suppression = -60):
        self.pulse = self.amplitude * signal.gausspulse(self.tspace - self.t_centre, self.frequency, self.bandwidth, suppression)
        self.spectrum = np.fft.fft(self.pulse)
        return self.pulse

    def get_spectrum(self): #NOTE: pulse must be defined before
        return self.spectrum