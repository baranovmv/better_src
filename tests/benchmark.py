import copy

from resamplers import *

import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann
from scipy.signal import stft
import matplotlib
import matplotlib.pyplot as plt
import unittest
from resamplers import *

class TestSrcMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000
        self.n_channels = 1
        self.src = Src(self.fs_in, self.fs_out, 1., self.n_channels)
        self.speex = SpeexResampler(self.fs_in, self.fs_out, 1., self.n_channels)
        # self.src = src_open(SRC_PROFILE_DEFAULT, self.n_channels-1,  self.fs_in, self.fs_out)
        error = ctypes.c_int()
        self.speex_state = speex_resampler_init(1, self.fs_in, self.fs_out, 10, error)
        self.initial_out_countdown = speex_resampler_get_output_latency(self.speex_state)
        self.initial_in_latency = speex_resampler_get_input_latency(self.speex_state)

        self.in_frame_size = self.in_frame_pos = min(self.initial_in_latency * self.n_channels, self.FrameSz);


        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        self.src.tear_down()
        self.speex.tear_down()
        # src_close(self.src)
        # speex_resampler_destroy(self.speex_state)

    def do_resample_src(self, sig_in, coeff=1.0, nchannels=1):
        times = self.src.do_resample(sig_in, coeff)
        times = np.array(times) / self.src.FrameSz * self.fs_in
        self.sig_out = self.src.sig_out
        self.sig_out_t = self.src.sig_out_t
        return np.mean(times).item()

    def do_resample_speex(self, sig_in, coeff=1.0):
        times = self.speex.do_resample(sig_in, coeff)
        times = np.array(times) / self.src.FrameSz * self.fs_in
        self.sig_out = self.speex.sig_out
        self.sig_out_t = self.speex.sig_out_t
        return np.mean(times).item()

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
        x = np.arange(-self.FrameSz*5, self.FrameSz*50)
        # x = np.ones(256) - 41
        t = np.arange(0, x.shape[0])
        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        # plt.plot(t, x, '.')
        self.do_resample_src(x)
        y_ref_src = finterp(self.sig_out_t)
        plt.plot(self.sig_out_t, y_ref_src-self.sig_out, '+', label="SRC")
        self.do_resample_speex(x)
        y_ref_speex = finterp(self.sig_out_t)
        plt.plot(self.sig_out_t, y_ref_speex-self.sig_out, 'x', label="SpeexDSP")

        plt.grid(True)
        plt.legend()
        plt.show()

        self.assertTrue(np.all(np.abs(y_ref_src-self.sig_out)/self.sig_out < 1e-3))


    def test_sinewave_upsample(self):
        fs = np.pi / 63
        n = np.arange(0, 1008*100)
        s_inp = np.sin(fs * n)
        finterp = lambda t: np.sin(fs*t)
        coef = 63/64
        src_time = self.do_resample_src(s_inp, coef)
        print(f"Src time: {src_time} ms")
        src_out = copy.deepcopy(self.sig_out)
        src_out_t = copy.deepcopy(self.sig_out_t)

        y_ref_src = finterp(src_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref_src)/self.sig_out, 0)
        inp_f, inp_db = self.avg_spect(s_inp, window_size=1008)
        src_f, src_db = self.avg_spect(src_out, window_size=1024)
        ref_f, ref_db = self.avg_spect(y_ref_src, window_size=1024)
        ref_max_db = np.max(ref_db)

        speex_time = self.do_resample_speex(s_inp, coef)
        print(f"Speex time: {speex_time} ms")
        self.sig_out = self.sig_out[4096:]
        self.sig_out_t = self.sig_out_t[4096:]
        speex_out = copy.deepcopy(self.sig_out)
        speex_out_t = copy.deepcopy(self.sig_out_t)

        y_ref_speex = finterp(self.sig_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref_speex)/self.sig_out, 0)
        speex_f, speex_db = self.avg_spect(self.sig_out, window_size=1024)

        # plt.subplot(211)
        # f, S = spectrum.compute_average_spectrum(y_ref_src-self.sig_out, window_size=1008, overlap=0)
        # plt.plot(f, S-ref_max_db, label="Input")
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()

        # plt.subplot(212)
        # plt.plot(*spectrum.compute_average_spectrum(s_inp, window_size=1024, overlap=0), label="Input")
        # plt.plot(*spectrum.compute_average_spectrum(y_ref_src, fs=48000/coef, window_size=1008, overlap=0), label="Ref SRC")
        # plt.plot(*spectrum.compute_average_spectrum(self.sig_out, fs=48000/coef, window_size=1008, overlap=0), label="Ref SRC")
        # plt.plot(*spectrum.compute_average_spectrum(speex_out, fs=48000/coef, window_size=1008, overlap=0), label="Ref SRC")
        # plt.plot(*spectrum.compute_average_spectrum(src_out, fs=48000/coef, window_size=1024, overlap=0), label="SRC")
        # plt.plot(*spectrum.compute_average_spectrum(speex_out, fs=48000/coef, window_size=1024, overlap=0), label="Speex")

        # plt.plot(inp_f, inp_db, label="Input")
        plt.plot(src_f, src_db, 'x', label="src")
        plt.plot(speex_f, speex_db, '+', label="speex")
        plt.plot(ref_f, ref_db, label="ref")
        plt.plot(inp_f, inp_db, label="input")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.legend()
        plt.show()

    def test_wgn_plot(self):
        n = np.arange(0, 1024*30)
        s_inp = np.random.randn(1024*30)
        # s_inp[:512] = 0
        # s_inp[195] = 1
        # t_ = 36.562500
        # dt_ = 1.125
        # sinc_step_ = 0.888888896
        coef = 0.75
        resamplers = {"SRC": self.src, "Speex": self.speex}
        for name, res in resamplers.items():
            _ = res.do_resample(s_inp, coef)
            freq, db = self.avg_spect(res.sig_out[4096:], window_size=1024)
            plt.plot(freq, db, label=name)

        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.show()


        # print(f"Time spent: {time_spent} ms")

        # x = self.sig_out[:512]

        # plt.subplot(211)
        # plt.plot(n, s_inp)
        # plt.plot(self.sig_out_t, self.sig_out)
        #
        # plt.subplot(212)
        # plt.plot(*spectrum.compute_average_spectrum(s_inp), label="Input")
        # plt.plot(*spectrum.compute_average_spectrum(self.sig_out, fs=48000/coef), label="Output")
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        # plt.show()


if __name__ == '__main__':
    matplotlib.use('TkAgg')
    unittest.main()
