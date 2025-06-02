#include "include/src.h"

#include <memory>

#include "SrcFactory.h"

src_t * src_open(SrcProfile_t profile, SrcNumChannels_t nchannels, const unsigned int in_fs, const unsigned int out_fs)
{
    Quality q;
    switch (profile)
    {
        case SRC_PROFILE_POOR:
            q = LOW;
            break;
        case SRC_PROFILE_MEDIUM:
        default:
            q = MEDIUM;
            break;
        case SRC_PROFILE_GOOD:
            q = HIGH;
            break;
    }
    ISrc* imp = SrcFactory::create(nchannels + 1, q, in_fs, out_fs);

    return imp->valid() ? (src_t*)imp : nullptr;
}

void src_set_scale(src_t * src, const float coeff)
{
    ISrc * imp = (ISrc*)src;
    imp->set_scaling(0, 0, coeff);
}

int src_push_samples(src_t * src, const float *samples, const unsigned int nsamples)
{
    ISrc * imp = (ISrc*)src;

    return imp->push(samples, nsamples) ? 1 : 0;
}

unsigned int src_pop_samples(src_t *src, float *out, const unsigned int max_out)
{
    ISrc * imp = (ISrc*)src;

    return imp->resample(out, max_out);
}

float src_left_to_process(const src_t * src)
{
    const ISrc * imp = (const ISrc*)src;
    return imp->left_2_process();
}

void src_close(src_t * src)
{
    ISrc * imp = (ISrc*)src;

    delete imp;
}

int src_push_samples(src_t *src, float *samples, unsigned int nsamples) {
    ISrc * imp = (ISrc*)src;

    return imp->push(samples, nsamples);
}


