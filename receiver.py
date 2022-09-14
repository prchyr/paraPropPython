# paraPropPython
# c. sbrocco, s. prohira
# A. Kyriacou

import util
import math as m
import numpy as np
from inspect import signature
from scipy.interpolate import interp1d
from scipy import signal

class receiver:
    """
    Parameters
    ----------
    x : float
        x position (m)
    z : float
        z position (m)
    """

    def __init__(self, x, z):
        self.x = x
        self.z = z

    def setup(self, freq, dt):
        """
        further setup of receiver using simulation parameters

        Parameters
        ----------
        freq : float array
            frequencies (GHz)
        dt : float
            time step (ns)
        """
        self.freq = freq
        self.spectrum = np.zeros(len(freq), dtype='complex')
        self.spectrum_plus = np.zeros(len(freq), dtype='complex')
        self.spectrum_minus = np.zeros(len(freq), dtype='complex')
        self.time = np.arange(0, dt * len(freq), dt)

    def add_spectrum_component(self, f, A):
        """
        adds the contribution of a frequency to the received signal spectrum

        Parameters
        ----------
        f : float
            corresponding frequencie (GHz)
        A : complex float
            complex amplitude of received siganl (V/m???)
        """
        i = util.findNearest(self.freq, f)
        self.spectrum[i] = A

    def add_spectrum_component_minus(self, f, A):
        """
                adds the contribution of a frequency to the received signal spectrum (for the reflected 'minus' signal

                Parameters
                ----------
                f : float
                    corresponding frequencie (GHz)
                A : complex float
                    complex amplitude of received siganl (V/m???)
                """
        i = util.findNearest(self.freq, f)
        self.spectrum_minus[i] = A

    def add_spectrum_component_plus(self, f, A):
        """
                adds the contribution of a frequency to the received signal spectrum (for the reflected 'minus' signal

                Parameters
                ----------
                f : float
                    corresponding frequencie (GHz)
                A : complex float
                    complex amplitude of received siganl (V/m???)
                """
        i = util.findNearest(self.freq, f)
        self.spectrum_plus[i] = A

    def get_spectrum(self):
        """
        gets received signal spectrum

        Returns
        -------
        1-d comlplex float array
        """
        return self.spectrum[:int(len(self.freq) / 2)]

    def get_signal(self):
        """
        gets received signal

        Returns
        -------
        1-d comlplex float array
        """
        return np.flip(util.doIFFT(self.spectrum))

    def get_signal_plus(self):
        '''
        get forward going (positive) signal

        Returns
        --------
        1-d complex float array
        '''
        return np.flip(util.doIFFT(self.spectrum_plus))

    def get_signal_minus(self):
        '''
         get backward going (minus) signal

        Returns
        --------
        1-d complex float array
        '''
        return np.flip(util.doIFFT(self.spectrum_minus))

    def get_frequency(self):
        """
        gets frequency array

        Returns
        -------
        1-d float array
        """
        return abs(self.freq)[:int(len(self.freq) / 2)]

    def get_time(self):
        """
        gets time array

        Returns
        -------
        1-d float array
        """
        return self.time