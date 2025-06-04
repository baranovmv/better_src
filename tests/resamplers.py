from src import *
from speexdsp_resampler import *

import ctypes
import numpy as np
from math import floor
import time

class SpeexResampler:
    FrameSz = 64
    Frame16Type = ctypes.c_float * FrameSz

    def __init__(self, fs_in_, fs_out_, coeff_, nchannels_):
        self.fs_in = fs_in_
        self.fs_out = fs_out_
        self.coeff = coeff_
        self.nchannels = nchannels_

        error = ctypes.c_int()
        self.speex_state = speex_resampler_init(self.nchannels, self.fs_in, self.fs_out, 4, error)

        self.initial_out_countdown = speex_resampler_get_output_latency(self.speex_state)
        self.initial_in_latency = speex_resampler_get_input_latency(self.speex_state)

        self.in_frame_size = self.in_frame_pos = min(self.initial_in_latency * self.nchannels, self.FrameSz)

        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])


    def tear_down(self):
        speex_resampler_destroy(self.speex_state)

    def set_scaling(self, coeff=1.0):
        max_numerator = 60000 # selected empirically
        base_frac = 10        # no more than 1 digit in fractional part

        base = round(max_numerator / max(self.fs_in, self.fs_out) * base_frac) / base_frac \
            if self.fs_in < max_numerator and self.fs_out < max_numerator \
            else 1.0

        ratio_num = round(self.fs_in * coeff * base)
        ratio_den = round(self.fs_out * base)

        err = speex_resampler_set_rate_frac(self.speex_state, ratio_num, ratio_den, round(self.fs_in * coeff), self.fs_out)

        self.speex_in_latency_diff = speex_resampler_get_input_latency(self.speex_state) - self.initial_in_latency

    def do_resample(self, sig_in, coeff=1.0):
        output_frame = self.Frame16Type()

        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

        self.set_scaling(coeff)
        dt = self.fs_in / self.fs_out * coeff
        t = 0
        time_spent_list = []

        idx = 0
        while idx < sig_in.shape[0]:
            x = sig_in[idx:idx+self.FrameSz].tolist()
            frame = self.Frame16Type(*x)

            in_len = ctypes.c_uint(len(x) // self.nchannels)
            out_len = ctypes.c_uint(self.FrameSz // self.nchannels)

            start_ts = time.time()
            err = speex_resampler_process_interleaved_float(self.speex_state, frame, in_len, output_frame, out_len)
            time_spent = time.time() - start_ts
            navailable = out_len.value * self.nchannels
            if navailable > 0:
                    time_spent_list.append(time_spent/navailable)

            # in_len.value now contains samples actually consumed (per-channel)
            idx += in_len.value * self.nchannels

            if self.initial_out_countdown > 0:
                n_samples = min(self.initial_out_countdown, out_len.value)
                navailable -= n_samples * self.nchannels
                self.initial_out_countdown -= n_samples

            if navailable == 0:
                continue

            sig_frame = np.array([output_frame[i] for i in range(navailable)])
            sig_frame_t = np.arange(0,dt * (navailable // self.nchannels), dt) + t
            t += dt * (navailable // self.nchannels)
            self.sig_out = np.concat((self.sig_out, sig_frame))
            sig_frame_t = np.repeat(sig_frame_t, self.nchannels)
            self.sig_out_t = np.concat((self.sig_out_t, sig_frame_t,))

        return time_spent_list


class Src:
    FrameSz = 64
    Frame16Type = ctypes.c_float * FrameSz

    def __init__(self, fs_in_, fs_out_, coeff_, nchannels_):
        self.fs_in = fs_in_
        self.fs_out = fs_out_
        self.coeff = coeff_
        self.nchannels = nchannels_

        self.src = src_open(SRC_PROFILE_MEDIUM, self.nchannels-1,  self.fs_in, self.fs_out)
        assert(self.src is not None)

        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

    def tear_down(self):
        src_close(self.src)

    def do_resample(self, sig_in, coeff=1.0):
        input_frame = self.Frame16Type()

        N = self.FrameSz
        npad = self.FrameSz - sig_in.shape[0] % self.FrameSz
        sig_split = [np.pad(sig_in[i:i + N], (0, N - sig_in[i:i + N].shape[0])) for i in range(0, sig_in.shape[0], N)]

        assert(src_set_scale(self.src, coeff) > 0)
        dt = self.fs_in / self.fs_out * coeff
        time_spent_list = []

        for i, x in enumerate(sig_split):
            frame = self.Frame16Type(*(x.tolist()))
            result = src_push_samples(self.src, frame, self.FrameSz)
            self.pushed += self.FrameSz
            while True:
                start_ts = time.time()
                navailable = src_pop_samples(self.src, input_frame, self.FrameSz)
                time_spent = time.time() - start_ts
                if navailable > 0:
                    time_spent_list.append(time_spent/navailable)
                # timestamp of the last sample
                tgap = self.pushed // self.nchannels - src_left_to_process(self.src)
                if navailable == 0:
                    break
                t = tgap - dt * (navailable // self.nchannels)
                sig_frame = np.array([input_frame[i] for i in range(navailable)])
                sig_frame_t = np.arange(0,dt * (navailable // self.nchannels), dt) + t
                self.sig_out = np.concat((self.sig_out, sig_frame))
                sig_frame_t = np.repeat(sig_frame_t, self.nchannels)
                self.sig_out_t = np.concat((self.sig_out_t, sig_frame_t,))

        npad_converted = floor(self.fs_out / self.fs_in / coeff * (npad // self.nchannels)) * self.nchannels
        if npad_converted > 0:
            self.sig_out_t = self.sig_out_t[:-npad_converted]
            self.sig_out = self.sig_out[:-npad_converted]
        return time_spent_list
