from src import *

import ctypes
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import unittest

class TestSrcMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz
    winlen = 33

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000
        self.src = src_open(SRC_PROFILE_DEFAULT, MONO,  self.fs_in, self.fs_out)
        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

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
        t = 0 if self.sig_out_t.shape[0] == 0 else self.sig_out_t[-1]

        for i, x in enumerate(sig_split):
            frame = self.Frame16Type(*(x.tolist()))
            result = src_push_samples(self.src, frame, self.FrameSz)
            self.assertEqual(1, result)
            self.pushed += self.FrameSz
            while True:
                navailable = src_pop_samples(self.src, input_frame, 16)
                # timestamp of the last sample
                tgap = self.pushed - src_left_to_process(self.src)
                if navailable == 0:
                    break
                t = tgap - dt * navailable
                sig_frame = np.array([input_frame[i] for i in range(navailable)])
                self.sig_out = np.concat((self.sig_out, sig_frame))
                sig_frame_t = np.arange(0,dt * navailable, dt) + t
                self.sig_out_t = np.concat((self.sig_out_t, sig_frame_t,))
                t = self.sig_out_t[-1] + dt

    def test_linear(self):
        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        # x = np.ones(256) - 41
        t = np.arange(0, x.shape[0])
        self.do_resample(x)

        finterp = interpolate.InterpolatedUnivariateSpline(t, x, k=3)
        # finterp = interp1d(t, x, kind='linear', assume_sorted=True)
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
        fs = np.pi / 8
        n = np.arange(0, 1024*30)
        s_inp = np.sin(fs * n)
        self.do_resample(s_inp, 1.7)

        finterp = interpolate.InterpolatedUnivariateSpline(n, s_inp, k=3)
        finterp = lambda t: np.sin(fs*t)
        y_ref = finterp(self.sig_out_t)
        self.assertTrue(np.all(np.abs(y_ref-self.sig_out)/self.sig_out < 1e-3))
        plt.plot(n, s_inp)
        plt.plot(self.sig_out_t, self.sig_out, 'o-')
        plt.plot(self.sig_out_t, y_ref, '+')
        plt.grid(True)
        plt.show()

if __name__ == '__main__':
    unittest.main()
