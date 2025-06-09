import copy

from resamplers import *

import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann, kaiser
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

        self.pushed = 0
        self.sig_out = np.array([])
        self.sig_out_t = np.array([])

        # To help clion to attach before the execution begins.
        time.sleep(0.5)

    def tearDown(self):
        pass
        # self.src.tear_down()
        # self.speex.tear_down()
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

    def avg_spect(self, x, fs=48000, window = np.ones(1024)):
        x = np.array(x)
        window_size = window.shape[0]

        # Compute STFT
        f, t, Zxx = stft(x, fs=fs, window=window, nperseg=window_size,
                         noverlap=0, return_onesided=True, boundary=None,
                         padded=False)
        spectra = np.abs(Zxx)
        avg_spectrum = np.mean(spectra, axis=1)
        avg_spectrum_db = 20 * np.log10(avg_spectrum + 1e-20)  # Add small value to avoid log(0)

        return f, avg_spectrum_db

    def generate_sine(self, freq, phase, duration, sr, t = None):
        if t is None:
            t = np.arange(0, duration, 1.0/sr)
        return np.sin(2 * np.pi * freq * t + phase)

    def generate_sweep(start_freq, end_freq, duration, sr):
        t = np.arange(0, duration, 1.0/sr)
        # Logarithmic sweep for better frequency coverage
        k = np.exp(np.log(end_freq/start_freq) / duration)
        phi = 2 * np.pi * start_freq * ((k**t - 1) / np.log(k))
        return np.sin(phi)

    def generate_complex_signal(duration, sr):
        # Generate a complex signal with multiple harmonics and transients
        t = np.arange(0, duration, 1.0/sr)
        signal = np.sin(2 * np.pi * 440 * t)  # A4 note
        signal += 0.5 * np.sin(2 * np.pi * 880 * t)  # First harmonic
        signal += 0.25 * np.sin(2 * np.pi * 1320 * t)  # Second harmonic

        # Add some percussive elements (short transients)
        for i in range(10):
            idx = int(i * sr * duration / 10)
            if idx + 100 < len(signal):
                signal[idx:idx+100] += np.exp(-np.arange(100)/20) * np.random.rand(100)

        # Normalize
        return signal / np.max(np.abs(signal))

    def calculate_thd(self, resampled, sr, fundamental_freq):
        """Calculate Total Harmonic Distortion"""

        window_size = round(sr / fundamental_freq)
        if resampled.shape[0] > 8192:
            n = 4096 // window_size
        else:
            n = 2
        window_size = window_size * n
        window = np.ones(window_size)

        # Compute STFT
        f, t, Zxx = stft(resampled, fs=sr, window=window, nperseg=window_size,
                         noverlap=0, return_onesided=True, boundary=None,
                         padded=False)
        spectra = np.abs(Zxx)
        avg_spectrum = np.mean(spectra[:, 2:], axis=1)
        fundamental_freq_bin = np.argmax(np.abs(avg_spectrum))
        fundamental_power = avg_spectrum[fundamental_freq_bin] ** 2
        harmonic_power = 0
        harmonics_bins = np.arange(fundamental_freq_bin*2, avg_spectrum.shape[0], step=fundamental_freq_bin)
        harmonic_power = np.sum(np.abs(avg_spectrum[harmonics_bins]) ** 2)
        thd = np.sqrt(harmonic_power / fundamental_power) * 100

        bin_list = np.insert(harmonics_bins, 0, fundamental_freq_bin)
        return thd, f[bin_list], avg_spectrum[bin_list], f, avg_spectrum

    def calculate_snr(self, reference, resampled):
        """Calculate Signal-to-Noise Ratio in dB"""
        # Make sure signals are the same length for comparison
        min_len = min(len(reference), len(resampled))
        reference = reference[:min_len]
        resampled = resampled[:min_len]

        # Calculate signal power
        signal_power = np.mean(reference ** 2)

        # Calculate noise power
        noise_power = np.mean((reference - resampled) ** 2)

        # Avoid division by zero
        if noise_power == 0:
            return float('inf')

        # Calculate SNR in dB
        snr = 10 * np.log10(signal_power / noise_power)
        return snr

    def calculate_spectral_flatness(self, signal, sr):
        """Calculate spectral flatness (Wiener entropy) - measure of noisiness vs. tonality"""
        # Using librosa implementation
        return librosa.feature.spectral_flatness(y=signal, sr=sr)[0].mean()

    def plot_spectrogram(self, signal, sr, title):
        plt.figure(figsize=(10, 4))
        D = librosa.amplitude_to_db(np.abs(librosa.stft(signal)), ref=np.max)
        librosa.display.specshow(D, y_axis='log', x_axis='time', sr=sr)
        plt.colorbar(format='%+2.0f dB')
        plt.title(title)
        return plt.gcf()

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


    def disable_test_sinewave_upsample(self):
        fs = np.pi / 63
        n = np.arange(0, 1008*100)
        s_inp = np.sin(fs * n)
        finterp = lambda t: np.sin(fs*t)
        coef = 63/64
        src_time = self.do_resample_src(s_inp, coef)
        src_out = copy.deepcopy(self.sig_out)
        src_out_t = copy.deepcopy(self.sig_out_t)

        y_ref_src = finterp(src_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref_src)/self.sig_out, 0)
        inp_f, inp_db = self.avg_spect(s_inp, window=np.ones(1008))
        src_f, src_db = self.avg_spect(src_out, window=np.ones(1024))
        ref_f, ref_db = self.avg_spect(y_ref_src, window=np.ones(1024))
        ref_max_db = np.max(ref_db)

        speex_time = self.do_resample_speex(s_inp, coef)
        self.sig_out = self.sig_out[4096:]
        self.sig_out_t = self.sig_out_t[4096:]
        speex_out = copy.deepcopy(self.sig_out)
        speex_out_t = copy.deepcopy(self.sig_out_t)

        y_ref_speex = finterp(self.sig_out_t)
        rel_err = np.where(self.sig_out > 1e-3, (self.sig_out-y_ref_speex)/self.sig_out, 0)
        speex_f, speex_db = self.avg_spect(self.sig_out, window=np.ones(1024))

        plt.subplot(211)
        plt.plot(n, s_inp, label="Inp")
        plt.plot(src_out_t, y_ref_src, label="Reference")
        plt.plot(src_out_t, src_out, label="SRC")
        plt.plot(speex_out_t, speex_out)

        plt.subplot(212)
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
        resamplers = { "Speex": SpeexResampler(self.fs_in, self.fs_out, 1., self.n_channels),
                       "SRC": Src(self.fs_in, self.fs_out, 1., self.n_channels),}
        for i, (name, res) in enumerate(resamplers.items()):
            # plt.subplot(len(resamplers), 1, i+1)
            _ = res.do_resample(s_inp, coef)
            freq, db = self.avg_spect(res.sig_out[4096:], window=np.ones(1024))
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

    def test_sine_harmonics(self):
        """Test harmonic distortion on pure sine waves"""

        results = {}
        # Test different frequency sines at different target sample rates
        # test_freqs = [100, 441, 960, 4410, 9600]
        target_srs = [44100, 48000, 96000]
        test_freqs = [2000]
        # target_srs = [24000]
        # resamplers = {"Speex": SpeexResampler}
        # resamplers = {"SRC": Src}
        resamplers = { "Speex": SpeexResampler, "SRC": Src,}
        duration = 10

        resampler_result = {}
        for resampler_name, resampler_ctr in resamplers.items():
            for freq in test_freqs:
                freq_results = {}

                for in_sr in [8000]:
                    input_signal = self.generate_sine(float(freq), 0, duration, in_sr)

                    for out_sr in [24000]:
                        coeff = 1.
                        resampler = resampler_ctr(in_sr, out_sr, coeff, 1)

                        # Skip if output SR is too low for the frequency (Nyquist)
                        if freq > out_sr / 2.1:  # Adding some margin
                            continue

                        # Generate a sine wave

                        times = resampler.do_resample(input_signal, coeff)
                        times = np.array(times) / resampler.FrameSz * self.fs_in
                        resampled_signal = resampler.sig_out
                        t_resampled = resampler.sig_out_t / in_sr
                        # t_resampled = np.arange(t_resampled[0], t_resampled[-1], step=1/out_sr)
                        reference_signal = self.generate_sine(float(freq), 0, duration, out_sr, t_resampled)

                        # Match lengths for comparison
                        min_len = min(len(resampled_signal), len(reference_signal))
                        resampled_signal = resampled_signal[:min_len]
                        reference_signal = reference_signal[:min_len]

                        # Calculate THD
                        thd, _, _, f, sf = self.calculate_thd(resampled_signal, out_sr, freq)
                        thdref, _, _, fref, sfref = self.calculate_thd(reference_signal, out_sr, freq)
                        # plt.plot(t_resampled[1:], (np.diff(t_resampled) - 1/out_sr), '+-')
                        plt.plot(f, 20*np.log10(sf+1e-20), label=resampler_name)

                        snr = self.calculate_snr(reference_signal, resampled_signal)

                        freq_results[out_sr] = {
                            'thd': thd,
                            'snr': snr
                        }

                        print(f"{resampler_name}: sine {freq}Hz @{in_sr} / {out_sr}Hz: THD={thd:.6f}%, SNR={snr:.2f}dB")

                results[freq] = freq_results
            resampler_result[resampler_name] = results

        plt.legend()
        plt.grid(True)
        plt.show()
        return results



if __name__ == '__main__':
    matplotlib.use('TkAgg')
    unittest.main()
