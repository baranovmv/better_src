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

class TestSrcMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000
        self.src = src_open(SRC_PROFILE_DEFAULT, MONO,  self.fs_in, self.fs_out)
        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        src_close(self.src)

    def do_resample(self, sig_in, coeff=1.0, nchannels=1):
        input_frame = self.Frame16Type()
        navailable = src_pop_samples(self.src, input_frame, self.FrameSz)
        self.assertEqual(0, navailable)

        N = self.FrameSz
        sig_split = [np.pad(sig_in[i:i + N], (0, N - sig_in[i:i + N].shape[0])) for i in range(0, sig_in.shape[0], N)]

        src_set_scale(self.src, coeff)
        dt = self.fs_in / self.fs_out * coeff
        time_spent = 0

        for i, x in enumerate(sig_split):
            frame = self.Frame16Type(*(x.tolist()))
            result = src_push_samples(self.src, frame, self.FrameSz)
            self.assertEqual(1, result)
            self.pushed += self.FrameSz
            while True:
                start_ts = time.time()
                navailable = src_pop_samples(self.src, input_frame, self.FrameSz)
                time_spent += time.time() - start_ts
                # timestamp of the last sample
                tgap = self.pushed // nchannels - src_left_to_process(self.src)
                if navailable == 0:
                    break
                self.assertTrue(navailable % nchannels == 0)
                t = tgap - dt * (navailable // nchannels)
                sig_frame = np.array([input_frame[i] for i in range(navailable)])
                sig_frame_t = np.arange(0,dt * (navailable // nchannels), dt) + t
                self.assertTupleEqual(sig_frame.shape, sig_frame_t.shape)
                self.sig_out = np.concat((self.sig_out, sig_frame))
                sig_frame_t = np.repeat(sig_frame_t, nchannels)
                self.sig_out_t = np.concat((self.sig_out_t, sig_frame_t,))

        return time_spent

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
        time.sleep(0.5)
        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        # x = np.ones(256) - 41
        t = np.arange(0, x.shape[0])
        self.do_resample(x)

        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        y_ref = finterp(self.sig_out_t)
        self.assertTrue(np.all(np.abs(y_ref-self.sig_out)/self.sig_out < 1e-3))

    def test_sinewave(self):
        fs = np.pi / 8
        n = np.arange(0, 96)
        s_inp = np.sin(fs * n)
        self.do_resample(s_inp)

        finterp = interpolate.InterpolatedUnivariateSpline(n, s_inp, k=3)
        y_ref = finterp(self.sig_out_t)
        self.assertTrue(np.all(np.abs(y_ref-self.sig_out)/self.sig_out < 1e-3))

    def test_sinewave_upsample(self):
        fs = np.pi / 63
        n = np.arange(0, 1008*30)
        s_inp = np.sin(fs * n)
        coef = 63/64
        self.do_resample(s_inp, coef)

        finterp = lambda t: np.sin(fs*t)
        y_ref = finterp(self.sig_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref)/self.sig_out, 0)

        _, avg_err_db = self.avg_spect(y_ref-self.sig_out, window_size=1008)
        _, avg_ref_db = self.avg_spect(y_ref, window_size=1008)
        ref_max_db = np.max(avg_ref_db)
        self.assertTrue(np.all(avg_err_db - ref_max_db < -120), np.max(avg_err_db))
        self.assertTrue(np.all(np.abs(rel_err) < 1e-3), np.max(np.abs(rel_err)))

        # n = n[:1024]
        # x = self.sig_out[n]
        # y = np.concat((x, x))
        #
        # plt.subplot(211)
        # plt.plot(n, x)
        # plt.plot(y)
        # plt.xlim([1030,1050])
        # plt.ylim([0.875, 1.1])
        # plt.plot(np.diff(self.sig_out_t))
        # plt.plot(self.sig_out_t, rel_err)
        # plt.plot(self.sig_out_t, self.sig_out)
        # f, S = spectrum.compute_average_spectrum(y_ref-self.sig_out, window_size=1008, overlap=0)
        # plt.plot(f, S-ref_max_db, label="Input")
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        #
        # plt.subplot(212)
        # plt.plot(*spectrum.compute_average_spectrum(s_inp, window_size=1008, overlap=0), label="Input")
        # plt.plot(*spectrum.compute_average_spectrum(y_ref, fs=48000/coef, window_size=1024, overlap=0), label="Input")
        # # plt.plot(*spectrum.compute_average_spectrum(self.sig_out, fs=48000/coef, window_size=1024, overlap=0), label="Input")
        # # # plt.plot(*spectrum.spect(y_ref[:1024], Fs=48000/coef))
        # plt.grid(True, alpha=0.3)
        # plt.tight_layout()
        # plt.show()


    def test_sinewave_downsample(self):
        time.sleep(0.5)
        fs = np.pi / 8
        n = np.arange(0, 1024*30)
        s_inp = np.sin(fs * n)
        coef = 1.125
        self.do_resample(s_inp, coef)

        finterp = lambda t: np.sin(fs*t)
        y_ref = finterp(self.sig_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref)/self.sig_out, 0)

        self.assertTrue(np.all(np.abs(rel_err) < 1e-3), np.max(np.abs(rel_err)))

    def test_wgn_plot(self):
        n = np.arange(0, 1024*30)
        s_inp = np.random.randn(1024*30)
        # s_inp[:512] = 0
        # s_inp[195] = 1
        # t_ = 36.562500
        # dt_ = 1.125
        # sinc_step_ = 0.888888896
        coef = 0.9
        time_spent = self.do_resample(s_inp, coef) * 1000
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

class TestSrcStereoMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz
    NChannels = 2

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000
        self.src = src_open(SRC_PROFILE_DEFAULT, STEREO,  self.fs_in, self.fs_out)
        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        src_close(self.src)

    def do_resample(self, sig_in, coeff=1.0):
        input_frame = self.Frame16Type()
        navailable = src_pop_samples(self.src, input_frame, self.FrameSz)
        self.assertEqual(0, navailable)

        N = self.FrameSz
        sig_split = [np.pad(sig_in[i:i + N], (0, N - sig_in[i:i + N].shape[0])) for i in range(0, sig_in.shape[0], N)]

        src_set_scale(self.src, coeff)
        dt = self.fs_in / self.fs_out * coeff
        time_spent = 0

        for i, x in enumerate(sig_split):
            frame = self.Frame16Type(*(x.tolist()))
            result = src_push_samples(self.src, frame, self.FrameSz)
            self.assertEqual(1, result)
            self.pushed += self.FrameSz
            while True:
                start_ts = time.time()
                navailable = src_pop_samples(self.src, input_frame, self.FrameSz)
                time_spent += time.time() - start_ts
                # timestamp of the last sample
                tgap = self.pushed // self.NChannels - src_left_to_process(self.src)
                if navailable == 0:
                    break
                self.assertTrue(navailable % self.NChannels == 0)
                t = tgap - dt * (navailable // self.NChannels)
                sig_frame = np.array([input_frame[i] for i in range(navailable)])
                sig_frame_t = np.arange(0, dt * (navailable // self.NChannels), dt) + t
                self.assertListEqual([x // self.NChannels for x in sig_frame.shape], list(sig_frame_t.shape))
                self.sig_out = np.concat((self.sig_out, sig_frame))
                sig_frame_t = np.repeat(sig_frame_t, self.NChannels)
                self.sig_out_t = np.concat((self.sig_out_t, sig_frame_t,))

        return time_spent

    def test_linear(self):
        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        y = -x + 100
        n = x.shape[0]
        x_stereo = np.repeat(x, 2)
        x_stereo[1::2] = y

        t = np.arange(0, n)
        self.do_resample(x_stereo, 0.5)

        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        x_ref = finterp(self.sig_out_t[0::2])
        finterp = interpolate.InterpolatedUnivariateSpline(t, y, k=3)
        y_ref = finterp(self.sig_out_t[1::2])
        self.assertTrue(np.all(np.abs(x_ref-self.sig_out[0::2])/self.sig_out[0::2] < 1e-3))
        self.assertTrue(np.all(np.abs(y_ref-self.sig_out[1::2])/self.sig_out[1::2] < 1e-3))

        # plt.subplot(211)
        # plt.plot(x, 'o-')
        # plt.plot(self.sig_out_t[0::2], self.sig_out[0::2], '+')
        #
        # plt.subplot(212)
        # plt.plot(y, 'o-')
        # plt.plot(self.sig_out_t[1::2], self.sig_out[1::2], '+')
        # plt.show()


if __name__ == '__main__':
    matplotlib.use('TkAgg')
    unittest.main()
