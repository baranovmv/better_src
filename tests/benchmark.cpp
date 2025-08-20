#include "src.h"
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <memory>
#include <random>

template <size_t N> void fill_wgn(std::array<float, N>& arr) {
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
constexpr SrcProfile_t profile = SRC_PROFILE_MEDIUM;

int main() {
    std::unique_ptr signal = std::make_unique<std::array<float, in_signal_duration> >();
    std::unique_ptr out_signal = std::make_unique<std::array<float, out_signal_duration> >();
    std::vector<double> time(out_signal_duration);
    double overall_time = 0;
    fill_wgn(*signal);

    std::cout << "Benchmarking SRC:" << std::endl;
    std::cout << "Input sample rate: " << in_fs << " Hz" << std::endl;
    std::cout << "Output sample rate: " << out_fs << " Hz" << std::endl;
    std::cout << "Quality: " << src_profile_to_str(profile) << std::endl;
    std::cout << "Frame size: " << frame_size << " samples" << std::endl;
    std::cout << "Total input duration: " << in_signal_duration / in_fs << " seconds" << std::endl;
    std::cout << std::endl;

    src_t* src = src_open(profile, MONO, in_fs, out_fs);
    size_t out_i = 0;
    for (size_t i = 0; i < in_signal_duration && out_i < out_signal_duration; i += frame_size) {
        assert(src_push_samples(src, &(*signal)[i], frame_size));
        size_t navailable = 0;
        do {
            const auto leftover = std::min(frame_size, out_signal->size() - out_i);
            auto start_time = std::chrono::high_resolution_clock::now();
            navailable = src_pop_samples(src, &(*out_signal)[out_i], leftover);
            const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::high_resolution_clock::now() - start_time)
                                        .count();
            if (navailable > frame_size / 2) {
                const double rate = duration * 1e-9 / navailable * out_fs;
                time.push_back(rate);
                overall_time += duration;
            }
            out_i += navailable;
        } while (navailable > 0 && out_i < out_signal_duration);
    }

    src_close(src);
    return 0;
}
