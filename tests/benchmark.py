import copy

from resamplers import *

import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann
from scipy.signal import stft
import matplotlib
import matplotlib.pyplot as plt
import unittest
import spectrum

class TestSrcMethods(unittest.TestCase):
    FrameSz = 16
    Frame16Type = ctypes.c_float * FrameSz

    def setUp(self):
        self.fs_in = 48000
        self.fs_out = 48000
        self.n_channels = 1
        self.src = src_open(SRC_PROFILE_DEFAULT, self.n_channels-1,  self.fs_in, self.fs_out)
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
        src_close(self.src)
        speex_resampler_destroy(self.speex_state)

    def do_resample_src(self, sig_in, coeff=1.0, nchannels=1):
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

    def set_scaling_speex(self, coeff=1.0):
        max_numerator = 60000 # selected empirically
        base_frac = 10        # no more than 1 digit in fractional part

        base = round(max_numerator / max(self.fs_in, self.fs_out) * base_frac) / base_frac \
            if self.fs_in < max_numerator and self.fs_out < max_numerator \
            else 1.0

        ratio_num = round(self.fs_in * coeff * base)
        ratio_den = round(self.fs_out * base)

        err = speex_resampler_set_rate_frac(self.speex_state, ratio_num, ratio_den, round(self.fs_in * coeff), self.fs_out)
        self.assertEqual(err, RESAMPLER_ERR_SUCCESS)

        self.speex_in_latency_diff = speex_resampler_get_input_latency(self.speex_state) - self.initial_in_latency

    def do_resample_speex(self, sig_in, coeff=1.0):
        output_frame = self.Frame16Type()

        self.sig_out = np.array([])
        self.sig_out_t = np.array([])
        N = self.FrameSz
        residual = []

        self.set_scaling_speex(coeff)
        dt = self.fs_in / self.fs_out * coeff
        t = 0
        time_spent = 0

        idx = 0
        while idx < sig_in.shape[0]:
            x = residual[:self.FrameSz]
            n_add_from_sig = self.FrameSz - len(x)
            x.extend(sig_in[idx:idx+n_add_from_sig].tolist())
            frame = self.Frame16Type(*x)
            remaining_in = ctypes.c_uint(self.FrameSz // self.n_channels)
            out_len = ctypes.c_uint(self.FrameSz // self.n_channels)

            start_ts = time.time()
            err = speex_resampler_process_interleaved_float(self.speex_state, frame, remaining_in, output_frame, out_len)
            time_spent += time.time() - start_ts

            residual = x[remaining_in.value * self.n_channels:]
            idx += remaining_in.value * self.n_channels
            navailable = out_len.value * self.n_channels
            if self.initial_out_countdown > 0:
                n_samples = min(self.initial_out_countdown, out_len.value)
                navailable -= n_samples * self.n_channels
                self.initial_out_countdown -= n_samples
            self.assertEqual(err, RESAMPLER_ERR_SUCCESS)

            if navailable == 0:
                continue

            sig_frame = np.array([output_frame[i] for i in range(navailable)])
            sig_frame_t = np.arange(0,dt * (navailable // self.n_channels), dt) + t
            t += dt * (navailable // self.n_channels)
            self.sig_out = np.concat((self.sig_out, sig_frame))
            sig_frame_t = np.repeat(sig_frame_t, self.n_channels)
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
        src_f, src_db = self.avg_spect(self.sig_out, window_size=1024)
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
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.legend()
        plt.show()


if __name__ == '__main__':
    matplotlib.use('TkAgg')
    unittest.main()
