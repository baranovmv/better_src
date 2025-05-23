# pragma once

class ISrc
{
public:
    virtual ~ISrc() = default;

    virtual bool set_scaling(size_t input_sample_rate, size_t output_sample_rate, float multiplier = 1.f) = 0;
    virtual bool push(const void *in, const size_t in_n) = 0;
    virtual size_t resample(void *out, const size_t out_sz) = 0;
    virtual bool valid() const = 0;
    virtual float left_2_process() const = 0;
};
