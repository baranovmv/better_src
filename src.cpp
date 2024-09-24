#include "include/src.h"
#include "converter.h"

src_t * src_open(SrcProfile_t profile, SrcNumChannels_t nchannels, const unsigned int in_fs, const unsigned int out_fs)
{
    Src* imp = nullptr;
    switch (profile) {
        case SRC_PROFILE_DEFAULT:
            imp = new Src(nchannels + 1, 90, in_fs, out_fs);
            break;
        default:
            return nullptr;
    }

    return imp->valid() ? (src_t*)imp : nullptr;
}

void src_set_scale(src_t * src, const float coeff)
{
    Src * imp = (Src*)src;
    imp->set_scaling(0, 0, coeff);
}

int src_push_samples(src_t * src, const float *samples, const unsigned int nsamples)
{
    Src * imp = (Src*)src;

    return imp->push(samples, nsamples) ? 1 : 0;
}

unsigned int src_pop_samples(src_t *src, float *out, const unsigned int max_out)
{
    Src * imp = (Src*)src;

    return imp->resample(out, max_out);
}

float src_left_to_process(const src_t * src)
{
    const Src * imp = (const Src*)src;
    return imp->left_2_process();
}

void src_close(src_t * src)
{
    Src * imp = (Src*)src;

    delete imp;
}

int src_push_samples(src_t *src, float *samples, unsigned int nsamples) {
    Src * imp = (Src*)src;

    return imp->push(samples, nsamples);
}


