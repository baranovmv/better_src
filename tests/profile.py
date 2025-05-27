from resamplers import *

import numpy as np
from scipy import interpolate
from scipy.signal.windows import hann
from scipy.signal import stft
import matplotlib
import matplotlib.pyplot as plt

def do_test(resampler):
    length = 48000 * 10
    n = np.arange(0, length)
    s_inp = np.random.randn(length)
    times_spent = resampler.do_resample(s_inp, 1.03)
    times_spent = np.array(times_spent) * 48000 * 1000

    return {"avg": times_spent.mean().item(),
            "var": times_spent.std().item(),
            "times": times_spent}


if __name__ == '__main__':
    matplotlib.use('TkAgg')

    src = Src(48000, 48000, coeff_=1.03, nchannels_=1)
    src_res = do_test(src)
    speex = SpeexResampler(48000, 48000, coeff_=1.03, nchannels_=1)
    speex_res = do_test(speex)

    print("Src:\t\tavg: {:.3f}\tstd: {:.3f}\t{}".format(src_res["avg"], src_res["var"], src_res["times"].shape[0] ))
    print("Speex:\t\tavg: {:.3f}\tstd: {:.3f}\t{}".format(speex_res["avg"], speex_res["var"], speex_res["times"].shape[0] ))

    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

    axs[0].hist(src_res["times"], bins = 71, range=(0, 15))
    axs[0].grid(True)
    axs[1].hist(speex_res["times"], bins = 71, range=(0, 15))
    axs[1].grid(True)
    plt.show()

