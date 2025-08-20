#include "src.h"
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <random>
#include <speex/speex_resampler.h>
#include <vector>
#include <algorithm>

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

constexpr double coeff = 1.03;
constexpr size_t in_fs = 48000;
constexpr size_t out_fs = 48000;
constexpr size_t in_signal_duration = in_fs * 100;
constexpr size_t out_signal_duration = out_fs * 100;
constexpr size_t frame_size = 64;
constexpr SrcProfile_t profile = SRC_PROFILE_MEDIUM;
constexpr int speex_quality = 4;
constexpr int channels = 1; // Mono

struct BenchmarkResult {
    double avg;
    double std;
    double overall_time;
    double realtimeness;
};

using Signal = std::array<float, in_signal_duration>;
using OutSignal = std::array<float, out_signal_duration>;

BenchmarkResult compute_stats(const std::vector<double>& time, double overall_time, size_t out_samples) {
    if (time.empty()) {
        return {0, 0, overall_time, 0};
    }
    // Copy and sort times
    std::vector<double> sorted_time = time;
    std::sort(sorted_time.begin(), sorted_time.end());
    // Exclude 15% smallest and 15% largest
    size_t n = sorted_time.size();
    size_t lower = static_cast<size_t>(n * 0.15);
    size_t upper = n - static_cast<size_t>(n * 0.15);
    if (upper <= lower) {
        lower = 0;
        upper = n;
    }
    double avg = 0;
    double std = 0;
    size_t count = upper > lower ? upper - lower : 0;
    if (count > 0) {
        for (size_t i = lower; i < upper; ++i)
            avg += sorted_time[i];
        avg /= count;
        for (size_t i = lower; i < upper; ++i)
            std += (sorted_time[i] - avg) * (sorted_time[i] - avg);
        std = sqrt(std / count);
    }
    double realtimeness =
        (out_samples > 0 ? overall_time * 1e-9 * static_cast<double>(out_fs) / static_cast<double>(out_samples) : 0);
    return { avg, std, overall_time, realtimeness };
}

template <typename QUALITY_T> void print_benchmark_header(const std::string& name, QUALITY_T quality) {
    std::cout << "Benchmarking " << name << ":" << std::endl;
    std::cout << "Input sample rate: " << in_fs << " Hz" << std::endl;
    std::cout << "Output sample rate: " << out_fs << " Hz" << std::endl;
    std::cout << "Quality: " << quality << std::endl;
    std::cout << "Frame size: " << frame_size << " samples" << std::endl;
    std::cout << "Total input duration: " << in_signal_duration / in_fs << " seconds" << std::endl;
}

void print_benchmark_results(const std::string& name, const BenchmarkResult& res) {
    std::cout << "\nBenchmark Results for " << name << ":\n";
    std::cout << "==============================\n";
    std::cout << "Average ms in second: " << res.avg * 1000. << std::endl;
    std::cout << "STD ms in second: " << res.std * 1000. << std::endl;
    std::cout << "Overall time: " << res.overall_time * 1e-9 << std::endl;
    std::cout << "Realtimeness: " << res.realtimeness << std::endl;
}

BenchmarkResult benchmark_src(const Signal& signal, OutSignal& out_signal) {
    std::vector<double> time;
    double overall_time = 0;
    src_t* src = src_open(profile, MONO, in_fs, out_fs);
    // Set scale for fractional coeff (like src_set_scale in Python)
    assert(src_set_scale(src, coeff) > 0);
    size_t out_samples = 0;
    for (size_t i = 0; i < in_signal_duration && out_samples < out_signal_duration; i += frame_size) {
        assert(src_push_samples(src, &signal[i], frame_size));
        size_t navailable = 0;
        do {
            const auto leftover = std::min(frame_size, out_signal.size() - out_samples);
            auto start_time = std::chrono::high_resolution_clock::now();
            navailable = src_pop_samples(src, &out_signal[out_samples], leftover);
            const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::high_resolution_clock::now() - start_time)
                                        .count();
            if (navailable > frame_size / 2) {
                const double rate = duration * 1e-9 / navailable * out_fs;
                time.push_back(rate);
                overall_time += duration;
            }
            out_samples += navailable;
        } while (navailable > 0 && out_samples < out_signal_duration);
    }
    src_close(src);
    return compute_stats(time, overall_time, out_samples);
}

BenchmarkResult benchmark_speex(const Signal& signal, OutSignal& out_signal) {
    std::vector<double> time;
    double overall_time = 0;
    int err = 0;
    SpeexResamplerState* resampler = speex_resampler_init(channels, in_fs, out_fs, speex_quality, &err);
    if (err != RESAMPLER_ERR_SUCCESS) {
        std::cerr << "Failed to initialize SpeexDSP resampler: " << err << std::endl;
        return { 0, 0, 0, 0 };
    }

    // --- Calculate ratio_num and ratio_den as in Python ---
    const int max_numerator = 60000;
    const int base_frac = 10;
    double base = (in_fs < max_numerator && out_fs < max_numerator)
        ? (std::round(max_numerator / std::max(in_fs, out_fs) * base_frac) / static_cast<double>(base_frac))
        : 1.0;
    spx_uint32_t ratio_num = static_cast<spx_uint32_t>(std::round(in_fs * coeff * base));
    spx_uint32_t ratio_den = static_cast<spx_uint32_t>(std::round(out_fs * base));
    // ------------------------------------------------------

    err = speex_resampler_set_rate_frac(
        resampler,
        ratio_num,
        ratio_den,
        static_cast<spx_uint32_t>(std::round(in_fs * coeff)),
        out_fs
    );
    if (err != RESAMPLER_ERR_SUCCESS) {
        std::cerr << "Failed to set rate fraction: " << err << std::endl;
        speex_resampler_destroy(resampler);
        return { 0, 0, 0, 0 };
    }

    size_t out_samples = 0, in_i = 0;
    float* in_pos = const_cast<float*>(signal.data());
    float* out_pos = out_signal.data();
    while (in_i < in_signal_duration && out_samples < out_signal_duration) {
        spx_uint32_t in_len = std::min(frame_size, in_signal_duration - in_i);
        spx_uint32_t out_len = std::min(static_cast<spx_uint32_t>(frame_size),
                                        static_cast<spx_uint32_t>(out_signal_duration - out_samples));
        auto start_time = std::chrono::high_resolution_clock::now();
        err = speex_resampler_process_interleaved_float(resampler, in_pos, &in_len, out_pos, &out_len);
        const double duration =
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time)
                .count();
        const size_t navailable = out_len * channels;
        if (navailable > frame_size / 2) {
            const double rate = duration * 1e-9 / navailable * out_fs;
            time.push_back(rate);
            overall_time += duration;
        }
        if (err != RESAMPLER_ERR_SUCCESS) {
            std::cerr << "Resampling error: " << err << std::endl;
            break;
        }
        out_samples += out_len * channels;
        in_pos += in_len * channels;
        out_pos += out_len * channels;
        in_i += in_len * channels;
    }
    speex_resampler_destroy(resampler);
    return compute_stats(time, overall_time, out_samples);
}

int main() {
    auto signal = std::make_unique<Signal>();
    auto out_signal_src = std::make_unique<OutSignal>();
    auto out_signal_speex = std::make_unique<OutSignal>();
    fill_wgn(*signal);

    print_benchmark_header("SpeexDSP resampler", speex_quality);
    BenchmarkResult speex_result = benchmark_speex(*signal, *out_signal_speex);
    print_benchmark_results("SpeexDSP", speex_result);

    std::cout << std::endl;
    print_benchmark_header("SRC", src_profile_to_str(profile));
    BenchmarkResult src_result = benchmark_src(*signal, *out_signal_src);
    print_benchmark_results("SRC", src_result);

    return 0;
}
