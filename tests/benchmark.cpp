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
    double avg = 0;
    for (auto x : time)
        avg += x;
    if (!time.empty())
        avg /= time.size();
    double std = 0;
    for (auto x : time)
        std += (x - avg) * (x - avg);
    if (!time.empty())
        std = sqrt(std / time.size());
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
