#include <iostream>
#include <math.h>
#include <vector>

#include "fixedpoint.h"

typedef float sample_t;

class Src
{
public:
    Src(size_t n_channels, const size_t win_len, const size_t in_fs, const size_t out_fs)
    : valid_(true)
    , in_fs_(in_fs)
    , out_fs_(out_fs)
    , n_channels_(n_channels)
    , win_len_(win_len)
    , win_len_effective_(win_len_)
    , win_len_effective_half_(float(win_len_effective_/2 + 1))
    , win_len_max_(3 * win_len)
    , middle_i_(win_len_max_ * n_channels_)
    , sinc_center_i_(win_len * window_interp_ / 2)
    , t_(0)
    , dt_(0)
    , sinc_step_(0)
    , delay_line_(win_len_max_ * 2 * n_channels_)
    , sinc_table_(win_len_ * window_interp_ + 1)
    , delay_line_i_(0)
    , delay_line_processed_i_((size_t)floorf((float)middle_i_ - win_len_effective_half_))
    , accum_low_(n_channels_)
    , accum_high_(n_channels_)
    {
        if (n_channels_ == 0) {
            valid_ = false;
            return;
        }
        if (win_len % 2 != 1 || win_len < 5) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }

        set_scaling(in_fs, out_fs, 1.f);
        fill_sinc_();
    }

    bool set_scaling(size_t input_sample_rate,
                     size_t output_sample_rate,
                     float multiplier = 1.f)
    {
        if (input_sample_rate == 0){
            input_sample_rate = in_fs_;
        }
        if(output_sample_rate == 0) {
            output_sample_rate = out_fs_;
        }
        in_fs_ = input_sample_rate;
        out_fs_ = output_sample_rate;

        const float new_scaling = float(input_sample_rate) / output_sample_rate * multiplier;

        // Filter out obviously invalid values.win_len_effective_half_
        if (new_scaling <= 0) {
            // TODO: error explanation.
            return false;
        } else if (new_scaling > 3) {
            // TODO: error explanation.
            return false;
        }

        sinc_step_ = cutoff_freq_ / std::max(1.f, new_scaling);
        win_len_effective_ = (size_t)ceilf((float)win_len_ / sinc_step_);
        win_len_effective_half_ = (float)win_len_ / 2.f / sinc_step_;
        dt_ = new_scaling;

        if (win_len_effective_ > win_len_max_) {
            return false;
        }

        valid_ = true;
        return true;
    }

    bool push(const sample_t *in, const size_t in_n)
    {
        if (available() == middle_i_ - 1) {
            // TODO: error explanation.
            return false;
        }

        size_t new_delay_line_i;
        std::copy(in, in + in_n * n_channels_, delay_line_.begin() + delay_line_i_);
        if (middle_i_ >= delay_line_i_ + in_n * n_channels_) {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |□□□□■■■■■■■■□□□□□□□ □□□□■■■■■■■■□□□□□□□□□□□|
            new_delay_line_i = delay_line_i_ + in_n * n_channels_;
            std::copy(in, in + in_n * n_channels_, delay_line_.begin() + middle_i_ + delay_line_i_);
        } else {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |■■■■□□□□□□□□□□□■■■■ ■■■■□□□□□□□□□□□□□□□■■■■|
            new_delay_line_i = (delay_line_i_ + in_n * n_channels_) - middle_i_;
            std::copy(in, in + middle_i_ - delay_line_i_, delay_line_.begin() + middle_i_ + delay_line_i_);
            std::copy(in + middle_i_ - delay_line_i_, in + new_delay_line_i, delay_line_.begin());
        }

        delay_line_i_ =  new_delay_line_i;

        return true;
    }

    inline size_t available() const
    {
        return dist_(delay_line_processed_i_, delay_line_i_);
    }

    size_t resample(sample_t *out, const size_t out_sz)
    {
        if (available() / n_channels_ < win_len_effective_) {
            return 0;
        }
        if (out_sz % n_channels_ != 0) {
            // TODO: error explanation.
            return 0;
        }

        size_t out_i = 0;
        while(available() / n_channels_ >= win_len_effective_ && out_i < out_sz) {
            const float cur_t = t_ >= win_len_effective_half_ ? t_ : t_ + float(middle_i_);
            const float win_start = cur_t - win_len_effective_half_;
            const float win_end = cur_t + win_len_effective_half_;

            size_t in_idx = (size_t) ceilf(win_start);
            const size_t in_idx_end = (size_t) floorf(win_end);
            const float sinc_start = ((float) sinc_center_i_) - (cur_t - (float) in_idx) * window_interp_ * sinc_step_;
            const size_t sinc_idx_step = size_t(sinc_step_ * window_interp_);
            double sinc_start_integral_;
            const float sinc_idx_residual = (float) modf(sinc_start, &sinc_start_integral_);

            size_t sinc_idx_low = (size_t) floorf(sinc_start);
            size_t sinc_idx_high = sinc_idx_low + 1;
            std::fill(accum_high_.begin(), accum_high_.end(), 0.f);
            std::fill(accum_low_.begin(), accum_low_.end(), 0.f);
            for (size_t i = in_idx; i < in_idx_end; i++) {
                for (size_t nchan = 0; nchan < n_channels_; nchan++) {
                    const size_t in_idx_i = i * n_channels_ + nchan;
                    accum_low_[nchan] += delay_line_[in_idx_i] * sinc_table_[sinc_idx_low];
                    accum_high_[nchan] += delay_line_[in_idx_i] * sinc_table_[sinc_idx_high];
                }
                sinc_idx_low += sinc_idx_step;
                sinc_idx_high = sinc_idx_low + 1;
            }

            for (size_t nchan = 0; nchan < n_channels_; ++nchan) {
                out[out_i * n_channels_ + nchan] = (accum_high_[nchan] - accum_low_[nchan]) * sinc_idx_residual
                                                   + accum_low_[nchan];
            }
            out_i++;
            t_ = t_ + dt_;
            if (t_ >= middle_i_) {
                t_ -= middle_i_;
            }
            // Update delay_line_processed_i_ if distance between it and t_ is greater
            // than half of window.
            if (t_ >= win_len_effective_half_) {
                delay_line_processed_i_ = (size_t)floorf(t_ - win_len_effective_half_);
            } else {
                delay_line_processed_i_ = (size_t)floorf(t_ + float(middle_i_) - win_len_effective_half_);
            }
        }

        return out_i;
    }

    bool valid() const
    {
        return valid_;
    }

private:
    static constexpr size_t SINC_INTERP_NBITS = 4;
    using time_t = FixedPoint<int32_t, 24>;
    using sinc_t = FixedPoint<int32_t, SINC_INTERP_NBITS>;

    bool valid_;
    size_t in_fs_;
    size_t out_fs_;
    const size_t n_channels_;
    const size_t win_len_;
    size_t win_len_effective_;
    float win_len_effective_half_; // Approximateion of win_len_effective_ / 2.
    const size_t win_len_max_;
    const size_t middle_i_;
    static constexpr size_t window_interp_{1 << SINC_INTERP_NBITS};
    static constexpr float sinc_unity_{1.0 / (double)window_interp_};
    const size_t sinc_center_i_;
    static constexpr float cutoff_freq_{1.f}; // TODO: return to 0.9f

    // Position of current output sample in terms of input samples (increments by 1/scaling), varies in [0, win_len_).
    float t_;
    // Increment of t_ -- 1/scaling.
    float dt_;
    float sinc_step_;

    std::vector<sample_t> delay_line_;
    std::vector<sample_t> sinc_table_;
    size_t delay_line_i_;
    // (t_ - win_len_) -- the first sample in the delay line which is still needed.
    size_t delay_line_processed_i_;

    std::vector<sample_t> accum_low_;
    std::vector<sample_t> accum_high_;


    inline size_t dist_(const size_t from, const size_t to) const
    {
        return from <= to ? to - from : to + middle_i_ - from;
    }

    static constexpr sample_t calc_sinc_(const sample_t x)
    {
        return sinf(M_PI * x) / M_PI / x;
    }

    void fill_sinc_()
    {
        double sinc_t = sinc_unity_;
        std::vector<double> window_log(sinc_table_.size());

        sinc_table_[sinc_center_i_] = 1.f;
        for (ssize_t i = 0; i < win_len_ * window_interp_ / 2; ++i) {
            const double window = 0.54
                                  - 0.46 * std::cos(2 * M_PI
                                               * ((double)(i - 1) / 2.0 / (double)sinc_table_.size() + 0.5));
            window_log[i] = window;
            sinc_table_[sinc_center_i_ - i - 1] = sinc_table_[sinc_center_i_ + i + 1] =
                    (float)(std::sin(M_PI * sinc_t) / M_PI / sinc_t * window);
            sinc_t += sinc_unity_;
        }
    }
};
