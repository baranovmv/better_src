// #include <immintrin.h>

#include <algorithm>

#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>

#include "fixedpoint.h"
#include "SrcInterface.h"

enum Quality
{
    LOW,
    MEDIUM,
    HIGH,
};

template <size_t N_CHANNELS, Quality QUALITY, typename sample_t = float, typename accum_t = sample_t>
class Src final : public ISrc
{
public:


    Src(const size_t in_fs, const size_t out_fs)
    : valid_(true)
    , started_(false)
    , in_fs_(in_fs)
    , out_fs_(out_fs)
    , win_len_effective_(WIN_LEN)
    , win_len_effective_half_(float(win_len_effective_/2 + 1))
    , win_len_max_(3 * WIN_LEN)
    , middle_i_(win_len_max_ * N_CHANNELS)
    , sinc_center_i_(WIN_LEN * window_interp_ / 2)
    , t_(ts_t::FromInteger(0))
    , t_counter_(ts_t::FromInteger(0))
    , t_win_begin_(ts_t::FromInteger(0))
    , dt_(ts_t::FromInteger(0))
    , sinc_step_(0)
    , delay_line_(win_len_max_ * 2 * N_CHANNELS)
    , sinc_table_(4 * WIN_LEN * window_interp_ + 1)
    , delay_line_i_(0)
    , delay_line_processed_i_(0)
    {
        static_assert(N_CHANNELS > 0 && N_CHANNELS <= 128, "n_channels must be between 1 and 128");

        if (WIN_LEN >= (1 << WINLEN_BITS)) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }

        if (N_CHANNELS == 0) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }
        if (WIN_LEN < 5) {
            // TODO: error explanation.
            valid_ = false;
            return;
        }

        set_scaling(in_fs, out_fs, 1.f);
    }

    bool set_scaling(size_t input_sample_rate,
                     size_t output_sample_rate,
                     float multiplier = 1.f) override
    {
        if (input_sample_rate == 0){
            input_sample_rate = in_fs_;
        }
        if(output_sample_rate == 0) {
            output_sample_rate = out_fs_;
        }
        in_fs_ = input_sample_rate;
        out_fs_ = output_sample_rate;

        const float new_scaling = float(input_sample_rate) / float(output_sample_rate) * multiplier;

        // Filter out obviously invalid values.win_len_effective_half_
        if (new_scaling <= 0 || new_scaling > 3) {
            valid_ = false;
            // TODO: error explanation.
            return false;
        }
        dt_ = ts_t(new_scaling);
        sinc_step_ = CUTOFF_FREQ / std::max(1.f, new_scaling);
        valid_ = fill_sinc_(sinc_step_);
        if (!started_) {
            t_ = win_len_effective_half_;
            t_counter_ = t_;
            t_win_begin_ = ts_t::FromInteger(0u);
            delay_line_processed_i_ = 0;
        } else {
            t_win_begin_ = t_ - win_len_effective_half_;
            delay_line_processed_i_ = size_t(ceil(t_win_begin_)) * N_CHANNELS;
        }

        return valid_;
    }

    bool push(const void *in, const size_t in_n) override
    {
        auto *in_samples = (const sample_t *)in;
        if (available() == middle_i_ - 1) {
            // TODO: error explanation.
            return false;
        }
        if (in_n % N_CHANNELS != 0) {
            // TODO: error explanation.
            return false;
        }

        size_t new_delay_line_i;
        std::copy(in_samples, in_samples + in_n, delay_line_.begin() + delay_line_i_);
        if (middle_i_ >= delay_line_i_ + in_n) {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |□□□□■■■■■■■■□□□□□□□ □□□□■■■■■■■■□□□□□□□□□□□|
            new_delay_line_i = delay_line_i_ + in_n;
            std::copy_n(in_samples, in_n, delay_line_.begin() + middle_i_ + delay_line_i_);
        } else {
            //  0               middle                   end
            //  ↓                  ↓                      ↓
            // |■■■■□□□□□□□□□□□■■■■ ■■■■□□□□□□□□□□□□□□□■■■■|
            new_delay_line_i = (delay_line_i_ + in_n) - middle_i_;
            std::copy_n(in_samples, middle_i_ - delay_line_i_, delay_line_.begin() + middle_i_ + delay_line_i_);
            std::copy_n(in_samples + middle_i_ - delay_line_i_, in_n, delay_line_.begin());
        }

        delay_line_i_ =  new_delay_line_i;

        return true;
    }

    [[nodiscard]] size_t available() const
    {
        return dist_(delay_line_processed_i_, delay_line_i_);
    }

    size_t resample(void *out, const size_t out_sz) override
    {
        auto *out_samples = (sample_t *)out;
        if (!started_ && (delay_line_i_ / N_CHANNELS < win_len_effective_)) {
            return 0;
        } else {
            started_ = true;
        }
        if (available() < win_len_effective_ * N_CHANNELS) {
            return 0;
        }
        if (out_sz % N_CHANNELS != 0) {
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
        while (available() > win_len_effective_ * N_CHANNELS
               && out_i < out_sz) {
            const auto offset = ts_t::FromInteger(delay_line_processed_i_ / N_CHANNELS) - t_win_begin_;
            sinc_t sinc_t_offset = sinc_t::FromInnerval(offset.get());

#if 1
            do_mac(sinc_t_offset, &out_samples[out_i]);
            out_i += N_CHANNELS;
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

            t_win_begin_ += dt_;
            if (t_win_begin_ >= ts_t(float(middle_i_ / N_CHANNELS))) {
                t_win_begin_ -= ts_t(float(middle_i_ / N_CHANNELS));
            }
            delay_line_processed_i_ = t_win_begin_.ceil() * ts_t::FromInteger(N_CHANNELS);
            counter_++;
        }
        t_ += ts_t::FromInteger(out_i / N_CHANNELS) * dt_;
        t_counter_ += ts_t::FromInteger(out_i) * dt_;
        while (t_ >= ts_t::FromInteger(middle_i_ / N_CHANNELS)) {
            t_ -= ts_t::FromInteger(middle_i_ / N_CHANNELS);
        }

        return out_i;
    }

    [[nodiscard]] bool valid() const override
    {
        return valid_;
    }

    [[nodiscard]] float left_2_process() const override
    {
        return dist_<float>(t_, float (delay_line_i_ / N_CHANNELS));
    }

private:

    static constexpr size_t get_winlen_bits(Quality q)
    {
        switch (q)
        {
        case LOW:
            return 6;
        case MEDIUM:
            return 8;
        case HIGH:
        default:
            return 10;
        }
    }

    static constexpr size_t get_winlen(Quality q)
    {
        switch (q)
        {
        case LOW:
            return 15;
        case MEDIUM:
            return 65;
        case HIGH:
        default:
            return 255;
        }
    }

    static constexpr size_t get_sinc_interp_bits(Quality q)
    {
        switch (q)
        {
        case LOW:
        case MEDIUM:
            return 4;
        case HIGH:
        default:
            return 5;
        }
    }

    static constexpr float get_cutoff_freq(Quality q)
    {
        switch (q)
        {
        case LOW:
            return 0.85;
        case MEDIUM:
            return 0.94;
        case HIGH:
        default:
            return 0.975;
        }
    }

    //! How many bits is enough to fit winlen_.
    static constexpr size_t WINLEN_BITS = get_winlen_bits(QUALITY);
    static constexpr size_t WIN_LEN = get_winlen(QUALITY);
    static constexpr size_t SINC_INTERP_NBITS = get_sinc_interp_bits(QUALITY);
    using ts_t = FixedPoint<uint32_t, uint64_t, 32 - WINLEN_BITS>;
    using sinc_t = FixedPoint<uint32_t, uint64_t, 32 - WINLEN_BITS - SINC_INTERP_NBITS>;
    static constexpr float CUTOFF_FREQ = get_cutoff_freq(QUALITY);

    bool valid_;
    bool started_;
    size_t in_fs_;
    size_t out_fs_;
    size_t win_len_effective_;
    ts_t win_len_effective_half_; // Approximateion of win_len_effective_ / 2.
    const size_t win_len_max_;
    const size_t middle_i_;
    static constexpr size_t window_interp_{1 << SINC_INTERP_NBITS};
    static constexpr float sinc_unity_{1.f / (float)window_interp_};
    size_t sinc_center_i_;

    // Position of current output sample in terms of input samples (increments by 1/scaling), varies in [0, win_len_).
    ts_t t_;
    ts_t t_counter_;
    ts_t t_win_begin_;
    // Increment of t_ -- 1/scaling.
    ts_t dt_;
    float sinc_step_;

    std::vector<sample_t> delay_line_;

    std::vector<sample_t> sinc_table_;
    size_t delay_line_i_;
    // (t_ - win_len_) -- the first sample in the delay line which is still needed.
    size_t delay_line_processed_i_;

    size_t counter_ = 0;

    template<class T>
    inline T dist_(const T from, const T to) const
    {
        return from <= to ? to - from : to + static_cast<T>(middle_i_) - from;
    }

    inline static constexpr sample_t calc_sinc_(const float x)
    {
        const float res = std::abs(x) < 1e-7 ? 1.f : static_cast<float>(std::sin(M_PI * x) / M_PI) / x;
        return static_cast<sample_t>(res);
    }

    inline static constexpr sample_t hann_win_(const size_t idx, const size_t len)
    {
        const auto n = double(len);
        const auto x = double(idx);
        return static_cast<sample_t>(0.5 - 0.5 * std::cos(2 * M_PI * x / n));
    }

    bool fill_sinc_(const float sinc_step)
    {
        win_len_effective_half_ = ts_t(float(WIN_LEN) / 2.f / sinc_step);

        sinc_center_i_ = size_t(ceilf(float(win_len_effective_half_) * window_interp_));
        win_len_effective_half_ = ts_t(float(sinc_center_i_) / window_interp_);
        win_len_effective_ = sinc_center_i_ * 2 / window_interp_;

        if (win_len_effective_ * 2 > win_len_max_) {
            // TODO: error explanation
            return false;
        }
        if (sinc_table_.size() < (win_len_effective_ + 2) * window_interp_) {
            // TODO: error explanation
            return false;
        }
        float sinc_idx = sinc_unity_;

        const sample_t amplitude = static_cast<sample_t>(sinc_step);
        sinc_table_[sinc_center_i_] = amplitude;
        for (ssize_t i = 1; i < sinc_center_i_; ++i) {
            const sample_t sinc_val = calc_sinc_(sinc_idx * sinc_step) * amplitude;
            const sample_t window = hann_win_(sinc_center_i_ + i, sinc_center_i_ * 2);
            sinc_table_[sinc_center_i_ - i] = sinc_table_[sinc_center_i_ + i] =
                     sinc_val * window;
            sinc_idx += sinc_unity_;
        }
        for (ssize_t i = sinc_center_i_*2; i < sinc_center_i_*2 + window_interp_ ; ++i) {
            const sample_t window = hann_win_(i, sinc_center_i_ * 2);
            sinc_table_[i] =
                    calc_sinc_(sinc_idx * sinc_step) *  amplitude * window;
            sinc_idx += sinc_unity_;
        }
        std::fill(sinc_table_.begin() + sinc_center_i_*2 + window_interp_,
                  sinc_table_.begin() + sinc_center_i_*2 + window_interp_*2,
                  static_cast<sample_t>(0));

        return true;
    }

    static void cubic_coef(float frac, std::array<float, 4> &interp)
    {
       /* Compute interpolation coefficients. I'm not sure whether this corresponds to cubic interpolation
       but I know it's MMSE-optimal on a sinc */
       interp[0] =  -0.16667f*frac + 0.16667f*frac*frac*frac;
       interp[1] = frac + 0.5f*frac*frac - 0.5f*frac*frac*frac;
       /*interp[2] = 1.f - 0.5f*frac - frac*frac + 0.5f*frac*frac*frac;*/
       interp[3] = -0.33333f*frac + 0.5f*frac*frac - 0.16667f*frac*frac*frac;
       /* Just to make sure we don't have rounding problems */
       interp[2] = 1.f-interp[0]-interp[1]-interp[3];
    }

    template<int N>
    static void lagrange_coef(float frac, std::array<float, N+1> &h)
    {
        // Initialize all elements to 1.0
        h.fill(1.0);

        // Compute Lagrange interpolation coefficients
        for (int k = 0; k <= N; ++k) {
            for (int n = 0; n <= N; ++n) {
                if (n != k) {
                    h[n] *= (frac - k) / (n - k);
                }
            }
        }
    }


    void do_mac(const sinc_t sinc_t_offset, sample_t *result)
    {
        if constexpr (false && N_CHANNELS == 1) {
            #if false &&  defined(__AVX512F__)
                // AVX-512 optimization for N_CHANNELS == 1
                constexpr size_t VEC_SIZE = 16; // 16 floats per __m512
                accum_t accum_low = 0;
                accum_t accum_high = 0;

                auto sinc_idx = sinc_t_offset.floor();
                size_t idx = delay_line_processed_i_;
                size_t end = delay_line_processed_i_ + win_len_effective_;
                size_t vec_end = idx + ((end - idx) / VEC_SIZE) * VEC_SIZE;

                __m512 acc_low = _mm512_setzero_ps();
                __m512 acc_high = _mm512_setzero_ps();

                for (; idx + VEC_SIZE <= vec_end; idx += VEC_SIZE, sinc_idx += window_interp_ * VEC_SIZE) {
                    // Load delay_line_ and sinc_table_ values
                    __m512 delay = _mm512_loadu_ps(&delay_line_[idx]);
                    __m512 sinc_l, sinc_h;

                    // Gather sinc_table_ values for low and high
                    __m512i sinc_indices = _mm512_set_epi32(
                        sinc_idx + window_interp_ * 15,
                        sinc_idx + window_interp_ * 14,
                        sinc_idx + window_interp_ * 13,
                        sinc_idx + window_interp_ * 12,
                        sinc_idx + window_interp_ * 11,
                        sinc_idx + window_interp_ * 10,
                        sinc_idx + window_interp_ * 9,
                        sinc_idx + window_interp_ * 8,
                        sinc_idx + window_interp_ * 7,
                        sinc_idx + window_interp_ * 6,
                        sinc_idx + window_interp_ * 5,
                        sinc_idx + window_interp_ * 4,
                        sinc_idx + window_interp_ * 3,
                        sinc_idx + window_interp_ * 2,
                        sinc_idx + window_interp_ * 1,
                        sinc_idx + window_interp_ * 0
                    );
                    __m512i sinc_indices_h = _mm512_add_epi32(sinc_indices, _mm512_set1_epi32(1));

                    sinc_l = _mm512_i32gather_ps(sinc_indices, sinc_table_.data(), 4);
                    sinc_h = _mm512_i32gather_ps(sinc_indices_h, sinc_table_.data(), 4);

                    acc_low = _mm512_fmadd_ps(delay, sinc_l, acc_low);
                    acc_high = _mm512_fmadd_ps(delay, sinc_h, acc_high);
                }

                // Horizontal sum of acc_low and acc_high
                float tmp[16];
                _mm512_storeu_ps(tmp, acc_low);
                for (int i = 0; i < 16; ++i) accum_low += tmp[i];
                _mm512_storeu_ps(tmp, acc_high);
                for (int i = 0; i < 16; ++i) accum_high += tmp[i];

                // Handle remaining elements
                for (; idx <= end; ++idx, sinc_idx += window_interp_) {
                    accum_low  += delay_line_[idx] * sinc_table_[sinc_idx];
                    accum_high += delay_line_[idx] * sinc_table_[sinc_idx + 1];
                }

                *result = sinc_t_offset.fract_linear_interp(accum_low, accum_high);
            #else
                accum_t accum_low_odd = 0;
                accum_t accum_high_odd = 0;
                accum_t accum_low_even = 0;
                accum_t accum_high_even = 0;

                auto sinc_idx_odd = sinc_t_offset.floor();
                auto sinc_idx_even = sinc_idx_odd + window_interp_;
                auto idx = delay_line_processed_i_;
                for (; idx < delay_line_processed_i_ + win_len_effective_; idx += 2) {
                    assert(sinc_idx_odd <= sinc_center_i_ * 2 + window_interp_);
                    assert(sinc_idx_even <= sinc_center_i_ * 2 + window_interp_);
                    accum_low_odd  +=  delay_line_[idx]     * sinc_table_[sinc_idx_odd];
                    accum_high_odd +=  delay_line_[idx]     * sinc_table_[sinc_idx_odd + 1];
                    accum_low_even  += delay_line_[idx + 1] * sinc_table_[sinc_idx_odd + window_interp_];
                    accum_high_even += delay_line_[idx + 1] * sinc_table_[sinc_idx_odd + window_interp_ + 1];
                    sinc_idx_odd += window_interp_ * 2;
                    // sinc_idx_even += window_interp_;
                }
                // idx -= 1;
                // sinc_idx_odd -= window_interp_;
                for (; idx <= delay_line_processed_i_ + win_len_effective_; idx += 1) {
                    accum_low_odd  += delay_line_[idx] * sinc_table_[sinc_idx_odd];
                    accum_high_odd += delay_line_[idx] * sinc_table_[sinc_idx_odd + 1];
                    sinc_idx_odd += window_interp_;
                }

                accum_low_odd += accum_low_even;
                accum_high_odd += accum_high_even;

                *result = sinc_t_offset.fract_linear_interp(accum_low_odd, accum_high_odd);
            #endif
        // Cubic interpolation between accumulators
        } else if (true) {
            std::array<accum_t, N_CHANNELS> accum_0;
            std::array<accum_t, N_CHANNELS> accum_1;
            std::array<accum_t, N_CHANNELS> accum_2;
            std::array<accum_t, N_CHANNELS> accum_3;
            std::array<float, 4> coef;
            cubic_coef(sinc_t_offset.fract()+1.f, coef);
            // lagrange_coef<3>(sinc_t_offset.fract()+1.f, coef);

            auto sinc_idx = sinc_t_offset.floor();
            auto idx = delay_line_processed_i_;
            for (auto nchan = 0; nchan < N_CHANNELS; ++nchan) {
                accum_0[nchan] = delay_line_[idx + nchan] * (sinc_idx > 0 ? sinc_table_[sinc_idx-1] : 0.f);
                accum_1[nchan] = delay_line_[idx + nchan] * sinc_table_[sinc_idx];
                accum_2[nchan] = delay_line_[idx + nchan] * sinc_table_[sinc_idx + 1];
                accum_3[nchan] = delay_line_[idx + nchan] * sinc_table_[sinc_idx + 2];
            }
            idx += N_CHANNELS;
            sinc_idx += window_interp_ - 1;
            size_t i = 1;
            const size_t i_end = win_len_effective_ * N_CHANNELS + 1;
            // for (; idx <= delay_line_processed_i_ + win_len_effective_ * N_CHANNELS;
                   // idx += N_CHANNELS, sinc_idx += window_interp_) {
            for (; i < i_end; ++i){
                assert(sinc_idx <= sinc_center_i_ * 2 + window_interp_);
                for (auto nchan = 0; nchan < N_CHANNELS; nchan++) {
                    accum_0[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx];
                    accum_1[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 1];
                    accum_2[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 2];
                    accum_3[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 3];
                }
                idx += N_CHANNELS;
                sinc_idx += window_interp_;
            }

            for (size_t nchan = 0; nchan < N_CHANNELS; nchan++) {
                *result++ =   coef[0] * accum_0[nchan]
                            + coef[1] * accum_1[nchan]
                            + coef[2] * accum_2[nchan]
                            + coef[3] * accum_3[nchan];
            }
        // Cubic interpolation reserve
        } else if (false) {
            std::array<accum_t, N_CHANNELS> accum_low_;
            std::array<accum_t, N_CHANNELS> accum_high_;

            std::fill(accum_high_.begin(), accum_high_.end(), 0.f);
            std::fill(accum_low_.begin(), accum_low_.end(), 0.f);

            auto sinc_idx = sinc_t_offset.floor();
            for (auto idx = delay_line_processed_i_; idx <= delay_line_processed_i_ + win_len_effective_ * N_CHANNELS; idx += N_CHANNELS) {
                assert(sinc_idx <= sinc_center_i_ * 2 + window_interp_);
                // Catmull-Rom spline coefficients
                const float y0 = sinc_idx > 1 ? sinc_table_[sinc_idx - 1] : 0.f; // y-1
                const float y1 = sinc_table_[sinc_idx]; // y0 (start point)
                const float y2 = sinc_table_[sinc_idx+1]; // y1 (end point)
                const float y3 = sinc_idx < sinc_center_i_ * 2 + window_interp_? sinc_table_[sinc_idx +2] : 0.f; // y2

                // Catmull-Rom cubic polynomial coefficients
                const float a0 = -0.5f * y0 + 1.5f * y1 - 1.5f * y2 + 0.5f * y3;
                const float a1 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
                const float a2 = -0.5f * y0 + 0.5f * y2;
                const float a3 = y1;

                const float t =  sinc_t_offset.fract();
                const float h = a0 * t * t * t + a1 * t * t + a2 * t + a3;
                // const float h  = a2 * t + a3;

                for (auto nchan = 0; nchan < N_CHANNELS; nchan++) {
                    accum_low_[nchan]  += delay_line_[idx + nchan] * h;
                    // accum_low_[nchan]  += delay_line_[idx + nchan] * sinc_table_[sinc_idx];
                    // accum_high_[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 1];
                }
                sinc_idx += window_interp_;
            }

            for (size_t nchan = 0; nchan < N_CHANNELS; nchan++) {
                *result++ = accum_low_[nchan];
                // *result++ = sinc_t_offset.fract_linear_interp(accum_low_[nchan], accum_high_[nchan]);
            }
        // Linear interpolation between accumulators
        } else {
            std::array<accum_t, N_CHANNELS> accum_low_;
            std::array<accum_t, N_CHANNELS> accum_high_;

            std::fill(accum_high_.begin(), accum_high_.end(), 0.f);
            std::fill(accum_low_.begin(), accum_low_.end(), 0.f);

            auto sinc_idx = sinc_t_offset.floor();
            for (auto idx = delay_line_processed_i_; idx <= delay_line_processed_i_ + win_len_effective_ * N_CHANNELS; idx += N_CHANNELS) {
                assert(sinc_idx <= sinc_center_i_ * 2 + window_interp_);
                for (auto nchan = 0; nchan < N_CHANNELS; nchan++) {
                    accum_low_[nchan]  += delay_line_[idx + nchan] * sinc_table_[sinc_idx];
                    accum_high_[nchan] += delay_line_[idx + nchan] * sinc_table_[sinc_idx + 1];
                }
                sinc_idx += window_interp_;
            }

            for (size_t nchan = 0; nchan < N_CHANNELS; nchan++) {
                *result++ = sinc_t_offset.fract_linear_interp(accum_low_[nchan], accum_high_[nchan]);
            }
        }
    }
};
