from resamplers import *

import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann
from scipy.signal import stft
import matplotlib
import matplotlib.pyplot as plt

fs_in = 48000
fs_out = 48000
coeff = 1.03


def do_test(resampler):
    length = fs_in * 100
    n = np.arange(0, length)
    s_inp = np.random.randn(length)
    # Pre-allocate output arrays for efficiency
    max_out_len = int(length * (fs_out // fs_in) * coeff * 1.05)  # allow for upsampling
    out_buf = np.zeros(max_out_len, dtype=np.float32)
    out_time_buf = np.zeros(max_out_len, dtype=np.float32)
    times_spent = resampler.do_resample(s_inp, coeff, out_buf, out_time_buf)
    times_spent = np.array(times_spent) * fs_in * 1000

    return {
        "avg": times_spent.mean().item(),
        "var": times_spent.std().item(),
        "times": times_spent,
    }


if __name__ == "__main__":
    matplotlib.use("TkAgg")

    speex = SpeexResampler(fs_in, fs_out, coeff_=coeff, nchannels_=1)
    speex_res = do_test(speex)
    src = Src(fs_in, fs_out, coeff_=coeff, nchannels_=1)
    src_res = do_test(src)

    print(
        "Src:\t\tavg: {:.3f}\tstd: {:.3f}\t{}".format(
            src_res["avg"], src_res["var"], src_res["times"].shape[0]
        )
    )
    print(
        "Speex:\t\tavg: {:.3f}\tstd: {:.3f}\t{}".format(
            speex_res["avg"], speex_res["var"], speex_res["times"].shape[0]
        )
    )

    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

    axs[0].hist(src_res["times"], bins=71, range=(0, 15))
    axs[0].grid(True)
    axs[1].hist(speex_res["times"], bins=71, range=(0, 15))
    axs[1].grid(True)
    plt.show()
