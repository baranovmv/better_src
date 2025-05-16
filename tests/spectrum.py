import numpy as np
from scipy.fft import fft, fftfreq
from math import *
import matplotlib.pyplot as plt
from scipy import signal
from numpy.random import randn


def spect(s, Fs=44100):
    """
    Use as follows: plt.plot(*spect(s)); plt.show()
    """
    spec = fft(s, norm="forward")
    xf = fftfreq(spec.size, d=1/Fs)

    n = round(s.shape[0]/2)
    return xf[:n], 20*np.log10(np.abs(spec[0:n]))


def pspect(s, Fs=44100):
    plt.plot(*spect(s, Fs))
    plt.grid(True)
    plt.show()


def semispect(s, Fs=44100):
    plt.semilogx(*spect(s, Fs))
    plt.grid(True, "both")
    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import windows, stft

def compute_average_spectrum(x, fs=48000, window_size=1024, kaiser_beta=12, overlap=0.5):
    """
    Compute the average spectrum of a signal using Kaiser windowing and STFT.

    Parameters:
    -----------
    x : array_like
        Input signal
    fs : int, optional
        Sampling rate in Hz (default: 48000)
    window_size : int, optional
        Size of the window in samples (default: 1024)
    kaiser_beta : float, optional
        Shape parameter of the Kaiser window (default: 12)
    overlap : float, optional
        Overlap ratio between consecutive windows (default: 0.5)

    Returns:
    --------
    freq : ndarray
        Frequency bins in Hz
    avg_spectrum : ndarray
        Average magnitude spectrum in dB
    """
    # Make sure input is numpy array
    x = np.array(x)

    # Create Kaiser window
    window = windows.kaiser(window_size, kaiser_beta)

    # Calculate the number of samples to step between windows (hop length)
    hop_length = int(window_size * (1 - overlap))

    # Compute STFT
    f, t, Zxx = stft(x, fs=fs, window=window, nperseg=window_size,
                     noverlap=window_size-hop_length, return_onesided=True)

    # Compute magnitude spectrum (normalized by window energy for proper scaling)
    window_energy = np.sum(window**2)
    spectra = np.abs(Zxx) / np.sqrt(window_energy)

    # Compute average spectrum across time frames
    avg_spectrum = np.mean(spectra, axis=1)

    # Convert to dB (with reference to 1.0)
    avg_spectrum_db = 20 * np.log10(avg_spectrum + 1e-10)  # Add small value to avoid log(0)

    return f, avg_spectrum_db

def plot_average_spectrum(freq, avg_spectrum, title="Average Spectrum"):
    """
    Plot the average spectrum.

    Parameters:
    -----------
    freq : ndarray
        Frequency bins in Hz
    avg_spectrum : ndarray
        Average magnitude spectrum in dB
    title : str, optional
        Plot title (default: "Average Spectrum")
    """
    plt.figure(figsize=(10, 6))
    plt.plot(freq, avg_spectrum)
    plt.grid(True, alpha=0.3)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dB)')
    plt.title(title)
    plt.xlim(0, freq[-1])
    plt.tight_layout()
    plt.show()

