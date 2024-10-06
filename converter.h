#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#include "fixedpoint.h"

typedef float sample_t;

class Src
{
public:
    Src(size_t n_channels, const size_t win_len, const size_t in_fs, const size_t out_fs)
    : valid_(true)
    , started_(false)
    , in_fs_(in_fs)
    , out_fs_(out_fs)
    , n_channels_(n_channels)
    , win_len_(win_len)
    , win_len_effective_(win_len_)
    , win_len_effective_half_(float(win_len_effective_/2 + 1))
    , win_len_max_(3 * win_len)
    , middle_i_(win_len_max_ * n_channels_)
    , sinc_center_i_(win_len * window_interp_ / 2)
    , t_(uint32_t(0))
    , t_counter_(uint32_t(0))
    , t_win_begin_(uint32_t(0))
    , dt_(uint32_t(0))
    , sinc_step_(0)
    , delay_line_(win_len_max_ * 2 * n_channels_)
    , sinc_table_(4 * win_len_ * window_interp_ + 1)
    , delay_line_i_(0)
    , delay_line_processed_i_(0)
    , accum_low_(n_channels_)
    , accum_high_(n_channels_)
    {
        if (win_len >= (1 << WINLEN_BITS)) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }

        if (n_channels_ == 0) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }
        if (win_len < 5) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }

        set_scaling(in_fs, out_fs, 1.f);
    }

    bool set_scaling(size_t input_sample_rate,
                     size_t output_sample_rate,
                     sample_t multiplier = 1.f)
    {
        if (input_sample_rate == 0){
            input_sample_rate = in_fs_;
        }
        if(output_sample_rate == 0) {
            output_sample_rate = out_fs_;
        }
        in_fs_ = input_sample_rate;
        out_fs_ = output_sample_rate;

        const sample_t new_scaling = sample_t(input_sample_rate) / sample_t(output_sample_rate) * multiplier;

        // Filter out obviously invalid values.win_len_effective_half_
        if (new_scaling <= 0 || new_scaling > 3) {
            valid_ = false;
            // TODO: error explanation.
            return false;
        }
        dt_ = time_t(new_scaling);
        sinc_step_ = cutoff_freq_ / std::max((sample_t)1., new_scaling);
        valid_ = fill_sinc_(sinc_step_);
        if (!started_) {
            t_ = win_len_effective_half_;
            t_counter_ = t_;
            t_win_begin_ = time_t(0u);
            delay_line_processed_i_ = 0;
        } else {
            t_win_begin_ = t_ - win_len_effective_half_;
            delay_line_processed_i_ = size_t(ceil(t_win_begin_)) * n_channels_;
        }

        return valid_;
    }

    bool push(const sample_t *in, const size_t in_n)
    {
        if (available() == middle_i_ - 1) {
            // TODO: error explanation.
            return false;
        }
        if (in_n % n_channels_ != 0) {
            // TODO: error explanation.
            return 0;
        }

        size_t new_delay_line_i;
        std::copy(in, in + in_n, delay_line_.begin() + delay_line_i_);
        if (middle_i_ >= delay_line_i_ + in_n) {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |□□□□■■■■■■■■□□□□□□□ □□□□■■■■■■■■□□□□□□□□□□□|
            new_delay_line_i = delay_line_i_ + in_n;
            std::copy(in, in + in_n, delay_line_.begin() + middle_i_ + delay_line_i_);
        } else {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |■■■■□□□□□□□□□□□■■■■ ■■■■□□□□□□□□□□□□□□□■■■■|
            new_delay_line_i = (delay_line_i_ + in_n) - middle_i_;
            std::copy(in, in + middle_i_ - delay_line_i_, delay_line_.begin() + middle_i_ + delay_line_i_);
            std::copy(in + middle_i_ - delay_line_i_, in + in_n, delay_line_.begin());
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
        if (!started_ && (delay_line_i_ / n_channels_ < win_len_effective_)) {
            return 0;
        } else {
            started_ = true;
        }
        if (available() < win_len_effective_ * n_channels_) {
            return 0;
        }
        if (out_sz % n_channels_ != 0) {
            // TODO: error explanation.
            return 0;
        }

        // Signal
        // 0                                            middle_i_
        // |------------------------------------------------|------------------------------------------------|
        //         t_win_begin_         t_      t_win_begin_+win_effectvie_
        //              ↓               ↓                   ↓
        // |            □■□□□■□□□■□□□■□□□■□□□■□□□■□□□■□□□■□□□
        // | 0...1...2...3...4...5...6...7...8...9..10..11... ... ... ... ... ... ... ... ... |
        // |             ↑                               ↑                                ↑
        //        dl_ln_processed_i_       dl_ln_processed_i_ + win_effective          delay_line_i_
        size_t out_i = 0;
        while (available() > win_len_effective_ * n_channels_
               && out_i < out_sz) {
            sinc_t sinc_t_offset = sinc_t(time_t(uint32_t(delay_line_processed_i_)) - t_win_begin_);

            std::fill(accum_high_.begin(), accum_high_.end(), 0.f);
            std::fill(accum_low_.begin(), accum_low_.end(), 0.f);

            auto sinc_idx = sinc_t_offset.floor();
#if 1
            for (auto idx = delay_line_processed_i_;
                idx <= delay_line_processed_i_ + win_len_effective_ * n_channels_;
                idx += n_channels_) {
                assert(sinc_idx <= sinc_center_i_ * 2 + window_interp_);
                do_mac_(sinc_idx, idx);
                sinc_idx += window_interp_;
            }
#elif 1
            if (t_counter_ >= 83.){
                fout.flush();
            }
            const double t_curr = t_ < win_len_effective_half_ ? t_ + middle_i_ : t_;
            const double t_win_begin = ceil(t_curr - win_len_effective_half_);
            double t_win_begin_counter = ceil(t_counter_ - win_len_effective_half_);
            sample_t sinc_t = t_win_begin - t_curr;
            size_t i1 = 0;
            size_t i2 = 0;
            for (auto idx = size_t(t_win_begin);
                idx <= size_t(t_curr + win_len_effective_half_);
                idx += n_channels_) {
                if (idx < t_curr) {
                    i1++;
                } else {
                    i2++;
                }

                const sample_t sinc_val = calc_sinc_(sinc_t/sinc_step_) / sinc_step_;
                const sample_t win_val = hann_win_(idx-size_t(t_win_begin), size_t(win_len_effective_half_ * 2));
                const sample_t sinc_coef = sinc_val * win_val;
                const sample_t delay_line_ref = std::sin(M_PI/8.*t_win_begin_counter);
                if (true || std::abs(delay_line_ref-delay_line_[idx]) > 1e-6){
                    accum_low_[0] += delay_line_ref * sinc_coef;
                } else {
                    accum_low_[0] += delay_line_[idx] * sinc_coef;
                }
                sinc_t += 1.;
                t_win_begin_counter += 1.;
            }
            if (i1 != i2){
                sinc_t_offset = 1;
            }
            sinc_t_offset = 0;
#else
            auto idx_bgn = delay_line_processed_i_;
            auto idx_end = delay_line_processed_i_ + win_len_effective_ * n_channels_;
            auto sinc_idx_end = sinc_idx + window_interp_*(win_len_effective_);
            do {
                do_mac_(sinc_idx, idx_bgn);
                do_mac_(sinc_idx_end, idx_end);

                sinc_idx += window_interp_;
                sinc_idx_end -= window_interp_;
                idx_bgn += n_channels_;
                idx_end -= n_channels_;
            } while(idx_end > idx_bgn);
            if (idx_end == idx_bgn)
                do_mac_(sinc_idx, idx_bgn);
#endif

            for (size_t nchan = 0; nchan < n_channels_; nchan++) {
                out[out_i++] = sinc_t_offset.fract_linear_interp(accum_high_[nchan], accum_low_[nchan]);
            }

            t_win_begin_ += dt_;
            if (t_win_begin_ >= sample_t (middle_i_ / n_channels_)) {
                t_win_begin_ -= time_t(sample_t(middle_i_ / n_channels_));
            }
            delay_line_processed_i_ = t_win_begin_.ceil();
            counter_++;
        }
        t_ += time_t(uint32_t(out_i / n_channels_)) * dt_;
        t_counter_ += time_t(uint32_t(out_i)) * dt_;
        while (t_ >= time_t(uint32_t(middle_i_ / n_channels_))) {
            t_ -= time_t(uint32_t(middle_i_ / n_channels_));
        }

        return out_i;
    }

    bool valid() const
    {
        return valid_;
    }

    float left_2_process() const
    {
        return dist_<float>(t_, float (delay_line_i_));
    }

private:
    static constexpr size_t WINLEN_BITS = 8; //! How many bits is enough to fit winlen_.
    static constexpr size_t SINC_INTERP_NBITS = 4;
    using time_t = FixedPoint<uint32_t, uint64_t, 32 - WINLEN_BITS>;
    using sinc_t = FixedPoint<uint32_t, uint64_t, 32 - WINLEN_BITS - SINC_INTERP_NBITS>;

    bool valid_;
    bool started_;
    size_t in_fs_;
    size_t out_fs_;
    const size_t n_channels_;
    const size_t win_len_;
    size_t win_len_effective_;
    time_t win_len_effective_half_; // Approximateion of win_len_effective_ / 2.
    const size_t win_len_max_;
    const size_t middle_i_;
    static constexpr size_t window_interp_{1 << SINC_INTERP_NBITS};
    static constexpr sample_t sinc_unity_{1.f / (sample_t)window_interp_};
    size_t sinc_center_i_;
    static constexpr sample_t cutoff_freq_{1.0f}; // TODO: return to 0.9f

    // Position of current output sample in terms of input samples (increments by 1/scaling), varies in [0, win_len_).
    time_t t_;
    time_t t_counter_;
    time_t t_win_begin_;
    // Increment of t_ -- 1/scaling.
    time_t dt_;
    sample_t sinc_step_;

    std::vector<sample_t> delay_line_;

    std::vector<sample_t> sinc_table_;
    size_t delay_line_i_;
    // (t_ - win_len_) -- the first sample in the delay line which is still needed.
    size_t delay_line_processed_i_;

    std::vector<sample_t> accum_low_;
    std::vector<sample_t> accum_high_;

    size_t counter_ = 0;

    template<class T>
    inline T dist_(const T from, const T to) const
    {
        return from <= to ? to - from : to + (T)middle_i_ - from;
    }

    inline static constexpr sample_t calc_sinc_(const sample_t x)
    {
        return std::abs(x) < 1e-7 ? 1.f : static_cast<sample_t>(std::sin(M_PI * x) / M_PI) / x;
    }

    inline static constexpr sample_t hann_win_(const size_t idx, const size_t len)
    {
        const auto n = double(len);
        const auto x = double(idx);
        return sample_t (0.5 - 0.5 * std::cos(2 * M_PI * x / n));
    }

    bool fill_sinc_(const sample_t sinc_step)
    {
        win_len_effective_half_ = time_t(float(win_len_) / 2.f / sinc_step_);

        if (win_len_effective_ * 2 > win_len_max_) {
            // TODO: error explanation
            return false;
        }

        sinc_center_i_ = size_t(ceilf(float(win_len_effective_half_) * window_interp_));
        win_len_effective_half_ = time_t(float(sinc_center_i_) / window_interp_);
        win_len_effective_ = sinc_center_i_ * 2 / window_interp_;
        if (sinc_table_.size() < (win_len_effective_ + 2) * window_interp_) {
            // TODO: error explanation
            return false;
        }
        sample_t sinc_idx = sinc_unity_;

        sinc_table_[sinc_center_i_] = 1.f;
        for (ssize_t i = 1; i < sinc_center_i_; ++i) {
            const sample_t window = hann_win_(sinc_center_i_ + i, sinc_center_i_ * 2);
            sinc_table_[sinc_center_i_ - i] = sinc_table_[sinc_center_i_ + i] =
                    calc_sinc_(sinc_idx / sinc_step) * window;
            sinc_idx += sinc_unity_;
        }
        for (ssize_t i = sinc_center_i_*2; i < sinc_center_i_*2 + window_interp_ ; ++i) {
            const sample_t window = hann_win_(i, sinc_center_i_ * 2);
            sinc_table_[i] =
                    calc_sinc_(sinc_idx / sinc_step) * window;
            sinc_idx += sinc_unity_;
        }
        std::fill(sinc_table_.begin() + sinc_center_i_*2 + window_interp_,
                  sinc_table_.begin() + sinc_center_i_*2 + window_interp_*2, 0);

        return true;
    }

    inline void do_mac_(const size_t sinc_idx, const size_t idx) {
        for (auto nchan = 0; nchan < n_channels_; nchan++) {
            accum_low_[nchan]  += delay_line_[idx + nchan] * sinc_table_[sinc_idx];
            accum_high_[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 1];
        }
    }
};
