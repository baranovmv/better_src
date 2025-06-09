#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <random>
#include "src.h"

template<size_t N>
void fill_wgn(std::array<float, N>& arr) {
    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Normal distribution for white noise (mean=0, stddev=1)
    std::normal_distribution<float> dist(0.0f, 1.0f);

    // Fill array with white noise
    for (auto& element : arr) {
        element = dist(gen);
    }
}

constexpr size_t in_fs = 8000;
constexpr size_t out_fs = 24000;
constexpr size_t in_signal_duration = in_fs * 100;
constexpr size_t out_signal_duration = out_fs * 100;
constexpr size_t frame_size = 64;

int main()
{
    std::unique_ptr signal = std::make_unique<std::array<float, in_signal_duration>>();
    std::unique_ptr out_signal = std::make_unique<std::array<float, out_signal_duration>>();
    fill_wgn(*signal);

    src_t * src = src_open(SRC_PROFILE_MEDIUM, MONO, in_fs, out_fs);
    size_t out_i = 0;
    for (size_t i = 0; i < in_signal_duration && out_i < out_signal_duration; i += frame_size) {
        assert(src_push_samples(src, &(*signal)[i], frame_size));
        size_t navailable = 0;
        do
        {
            const auto leftover = std::min(frame_size, out_signal->size() - out_i);
            navailable = src_pop_samples(src, &(*out_signal)[out_i], leftover);
            out_i += navailable;
        } while (navailable > 0 && out_i < out_signal_duration);
    }

    src_close(src);
    return 0;
}
