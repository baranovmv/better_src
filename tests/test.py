from src import *

import ctypes
import numpy as np
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
                if navailable == 0:
                    break
                self.sig_out = np.concat((self.sig_out, np.array([input_frame[i] for i in range(navailable)])))
                self.sig_out_t = np.concat((self.sig_out_t, np.arange(t, t + dt * navailable, dt),))
                t = self.sig_out_t[-1] + dt

    def test_flow(self):
        x = np.arange(-self.FrameSz*2.5, self.FrameSz*2.5)
        t = np.arange(0, x.shape[0])
        self.do_resample(x)

        idx = np.isin(t, self.sig_out_t)
        t = t[idx]
        self.assertTrue(np.all(t == self.sig_out_t))
        self.assertTrue(np.all(np.abs(x[idx] - self.sig_out) < 1e-10))

        # plt.plot(t, x)
        # plt.plot(self.sig_out_t, self.sig_out)
        # plt.show()

if __name__ == '__main__':
    unittest.main()