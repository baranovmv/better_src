from src import *

import ctypes
import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann
from scipy.signal import stft
import matplotlib
import matplotlib.pyplot as plt
import time
import unittest
import spectrum

from resamplers import Src

class TestSrcMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        pass

    def avg_spect(self,x, fs=48000, window_size=1024):
        x = np.array(x)
        window = np.ones(window_size)

        # Compute STFT
        f, t, Zxx = stft(x, fs=fs, window=window, nperseg=window_size,
                         noverlap=0, return_onesided=True, boundary=None,
                         padded=False)
        spectra = np.abs(Zxx)
        avg_spectrum = np.mean(spectra, axis=1)
        avg_spectrum_db = 20 * np.log10(avg_spectrum + 1e-20)  # Add small value to avoid log(0)

        return f, avg_spectrum_db

    def test_linear(self):
        resampler = Src(48000, 48000, 1., 1)
        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        # x = np.ones(256) - 41
        t = np.arange(0, x.shape[0])
        _ = resampler.do_resample(x, 1.0)
        # self.do_resample(x)

        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        y_ref = finterp(resampler.sig_out_t)
        self.assertTrue(np.all(np.abs(y_ref-resampler.sig_out)/resampler.sig_out < 1e-3))

    def test_sinewave(self):
        freqs = [8000, 16000, 24000, 44100, 48000, 96000]
        sig_freq = [0.01, 0.1, 0.5, 0.7]
        for in_freq in [8000]:
            for out_freq in [24000]:
                for freq_coef in [0.5]:
                    resampler = Src(in_freq, out_freq, 1., 1)
                    fs = freq_coef * in_freq / 2
                    duration = 0.1
                    t = np.arange(0, duration, 1.0 / in_freq)
                    s_inp = np.sin(2*np.pi * fs * t)
                    resampler.do_resample(s_inp, 1.)

                    finterp = lambda t: np.sin(2*np.pi * fs * t)
                    y_ref = finterp(resampler.sig_out_t / in_freq)
                    self.assertTrue(np.all(np.abs(y_ref-resampler.sig_out)/resampler.sig_out < 1e-3), f"In sr: {in_freq}, Out sr: {out_freq}, freq_coef: {freq_coef}")

    def test_sinewave_upsample(self):
        resampler = Src(48000, 48000, 1., 1)
        fs = np.pi / 63
        n = np.arange(0, 1008*30)
        s_inp = np.sin(fs * n)
        coef = 63/64
        resampler.do_resample(s_inp, coef)

        finterp = lambda t: np.sin(fs*t)
        y_ref = finterp(resampler.sig_out_t)
        rel_err = np.where(resampler.sig_out > 1e-3, (resampler.sig_out-y_ref)/resampler.sig_out, 0)

        _, avg_err_db = self.avg_spect(y_ref-resampler.sig_out, window_size=1008)
        _, avg_ref_db = self.avg_spect(y_ref, window_size=1008)
        ref_max_db = np.max(avg_ref_db)
        self.assertTrue(np.all(avg_err_db - ref_max_db < -100), np.max(avg_err_db - ref_max_db))
        self.assertTrue(np.all(np.abs(rel_err) < 1e-3), np.max(np.abs(rel_err)))

        # n = n[:1024]
        # x = resampler.sig_out[n]
        # y = np.concat((x, x))
        #
        # plt.subplot(211)
        # plt.plot(n, x)
        # plt.plot(y)
        # plt.xlim([1030,1050])
        # plt.ylim([0.875, 1.1])
        # plt.plot(np.diff(resampler.sig_out_t))
        # plt.plot(resampler.sig_out_t, rel_err)
        # plt.plot(resampler.sig_out_t, resampler.sig_out)
        # f, S = spectrum.compute_average_spectrum(y_ref-resampler.sig_out, window_size=1008, overlap=0)
        # plt.plot(f, S-ref_max_db, label="Input")
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        #
        # plt.subplot(212)
        # plt.plot(*spectrum.compute_average_spectrum(s_inp, window_size=1008, overlap=0), label="Input")
        # plt.plot(*spectrum.compute_average_spectrum(y_ref, fs=48000/coef, window_size=1024, overlap=0), label="Input")
        # # plt.plot(*spectrum.compute_average_spectrum(resampler.sig_out, fs=48000/coef, window_size=1024, overlap=0), label="Input")
        # # # plt.plot(*spectrum.spect(y_ref[:1024], Fs=48000/coef))
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        # plt.show()


    def test_sinewave_downsample(self):
        time.sleep(0.5)
        resampler = Src(48000, 48000, 1., 1)
        fs = np.pi / 8
        n = np.arange(0, 1024*30)
        s_inp = np.sin(fs * n)
        coef = 1.125
        resampler.do_resample(s_inp, coef)

        finterp = lambda t: np.sin(fs*t)
        y_ref = finterp(resampler.sig_out_t)
        rel_err = np.where(resampler.sig_out > 1e-3, (resampler.sig_out-y_ref)/resampler.sig_out, 0)

        self.assertTrue(np.all(np.abs(rel_err) < 1e-3), np.max(np.abs(rel_err)))

    def test_wgn_plot(self):
        resampler = Src(48000, 48000, 1., 1)

        n = np.arange(0, 1024*30)
        s_inp = np.random.randn(1024*30)
        # s_inp[:512] = 0
        # s_inp[195] = 1
        # t_ = 36.562500
        # dt_ = 1.125
        # sinc_step_ = 0.888888896
        coef = 0.9
        time_spent = resampler.do_resample(s_inp, coef) * 1000
        # print(f"Time spent: {time_spent} ms")

        # x = resampler.sig_out[:512]

        # plt.subplot(211)
        # plt.plot(n, s_inp)
        # plt.plot(resampler.sig_out_t, resampler.sig_out)
        #
        # plt.subplot(212)
        # plt.plot(*spectrum.compute_average_spectrum(s_inp), label="Input")
        # plt.plot(*spectrum.compute_average_spectrum(resampler.sig_out, fs=48000/coef), label="Output")
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        # plt.show()

class TestSrcStereoMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz
    NChannels = 2

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        pass

    def test_linear(self):
        resampler = Src(48000, 48000, 1., 2)

        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        y = -x + 100
        n = x.shape[0]
        x_stereo = np.repeat(x, 2)
        x_stereo[1::2] = y

        t = np.arange(0, n)
        resampler.do_resample(x_stereo, 0.5)

        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        x_ref = finterp(resampler.sig_out_t[0::2])
        finterp = interpolate.InterpolatedUnivariateSpline(t, y, k=3)
        y_ref = finterp(resampler.sig_out_t[1::2])
        self.assertTrue(np.all(np.abs(x_ref-resampler.sig_out[0::2])/resampler.sig_out[0::2] < 1e-3))
        self.assertTrue(np.all(np.abs(y_ref-resampler.sig_out[1::2])/resampler.sig_out[1::2] < 1e-3))

        # plt.subplot(211)
        # plt.plot(x, 'o-')
        # plt.plot(resampler_t[0::2], resampler[0::2], '+')
        #
        # plt.subplot(212)
        # plt.plot(y, 'o-')
        # plt.plot(resampler_t[1::2], resampler[1::2], '+')
        # plt.show()


if __name__ == '__main__':
    matplotlib.use('TkAgg')
    unittest.main()
