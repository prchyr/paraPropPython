import util
import math as m
import numpy as np
from inspect import signature
from scipy.interpolate import interp1d
from scipy import signal
# A. Kyriacou

class tx_signal:
    def __init__(self, frequency, bandwidth, t_centre, dt, tmax, amplitude = 1, noise_amplitude = 0, freqMin=None, freqMax=None):
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
        self.noise_amplitude = noise_amplitude
        if freqMax == None:
            self.freqMax = self.frequency + self.bandwidth/2
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

    def do_impulse_response(self, IR, IR_freq):
        spectrum_shift = np.fft.fftshift(self.spectrum)
        freq_shift = np.fft.fftshift(self.freq_space)
        nFreq = len(freq_shift)
        IR_fmin = min(IR_freq)
        IR_fmax = max(IR_freq)
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
                k_pos = util.findNearest(IR_freq, freq_plus)
                self.IR_positive[i] = IR[k_pos]
                self.IR_negative[j_neg] = IR[k_pos]
        self.IR[nMid:] = self.IR_positive
        self.IR[:nMid] = np.flip(self.IR_negative)

        spectrum_shift *= self.IR
        self.spectrum = np.fft.ifftshift(spectrum_shift)
        self.pulse = np.fft.ifft(self.spectrum)

        return self.pulse

    def add_gaussian_noise(self):
        nSamples = len(self.spectrum)
        noise_amplitude = self.noise_amplitude
        if noise_amplitude > 0:
            self.noise = noise_amplitude * np.random.normal(0, noise_amplitude, nSamples)
            self.spectrum += self.noise
            self.pulse = np.fft.ifft(self.spectrum)

