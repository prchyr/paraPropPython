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

    def do_impulse_response(self, IR_freq, IR_data):
        '''
        Applies Convolution of received amplitude and antenna's impulse response function
        -> Simply applies
        '''
        self.IR = np.ones(self.nSamples)
        self.IR_freq = IR_freq
        self.IR_data = IR_data

        spectrum_shift = np.fft.fftshift(self.spectrum)
        freq_shift = np.fft.fftshift(self.freq_space)
        nFreq = len(freq_shift)
        IR_fmin = min(self.IR_freq)
        IR_fmax = max(self.IR_freq)
        nMid = util.findNearest(freq_shift, 0)

        freq_shift_positive = freq_shift[nMid:]
        freq_shift_negative = -1 * np.flip(freq_shift[:nMid])

        nFreq_pos = len(freq_shift_positive)
        nFreq_neg = len(freq_shift_positive)
        self.IR = np.zeros(nFreq)
        self.IR_positive = np.zeros(nFreq_pos)
        self.IR_negative = np.zeros(nFreq_neg)
        for i in range(nFreq_pos):
            freq_plus = freq_shift_positive[i]
            j_neg = util.findNearest(freq_shift_negative, freq_plus)
            freq_neg = freq_shift_negative[j_neg]
            if freq_plus < IR_fmin:
                self.IR_positive[i] = 0
                self.IR_negative[j_neg] = 0
            elif freq_plus > IR_fmax:
                self.IR_positive[i] = 0
                self.IR_negative[j_neg] = 0
            elif freq_plus <= IR_fmax and freq_plus >= IR_fmin:
                print('reset')
                k_pos = util.findNearest(self.IR_freq, freq_plus)
                self.IR_positive[i] = self.IR_data[k_pos]
                self.IR_negative[j_neg] = self.IR_data[k_pos]
        self.IR[nMid:] = self.IR_positive
        self.IR[:nMid] = np.flip(self.IR_negative)

        self.spectrum *= self.IR
        self.spectrum_plus *= self.IR
        self.spectrum_minus *= self.IR

    def get_impulse_response(self):
        return self.IR

    def add_gaussian_noise(self, noise_amplitude):
        self.noise_amplitude = noise_amplitude
        nSamples = len(self.spectrum)
        if self.noise_amplitude > 0:
            self.noise = self.noise_amplitude * np.random.normal(0, self.noise_amplitude, nSamples)
            self.spectrum += self.noise