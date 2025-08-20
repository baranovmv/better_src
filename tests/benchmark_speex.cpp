#include "src.h"
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <speex/speex_resampler.h>

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
constexpr size_t in_signal_duration = in_fs * 100;   // 100 seconds of audio
constexpr size_t out_signal_duration = out_fs * 100; // 100 seconds of output
constexpr size_t frame_size = 64;
constexpr int channels = 1; // Mono
constexpr int quality = 4;  // Medium quality (0-10, where 10 is best)

int main() {
    std::unique_ptr<std::array<float, in_signal_duration> > signal =
        std::make_unique<std::array<float, in_signal_duration> >();
    std::unique_ptr<std::array<float, out_signal_duration> > out_signal =
        std::make_unique<std::array<float, out_signal_duration> >();
    std::vector<double> time(out_signal_duration);

    fill_wgn(*signal);

    // Initialize SpeexDSP resampler
    int err = 0;
    SpeexResamplerState* resampler = speex_resampler_init(channels, in_fs, out_fs, quality, &err);

    if (err != RESAMPLER_ERR_SUCCESS) {
        std::cerr << "Failed to initialize SpeexDSP resampler: " << err << std::endl;
        return 1;
    }

    std::cout << "Benchmarking SpeexDSP resampler:" << std::endl;
    std::cout << "Input sample rate: " << in_fs << " Hz" << std::endl;
    std::cout << "Output sample rate: " << out_fs << " Hz" << std::endl;
    std::cout << "Quality: " << quality << std::endl;
    std::cout << "Frame size: " << frame_size << " samples" << std::endl;
    std::cout << "Total input duration: " << in_signal_duration / in_fs << " seconds" << std::endl;
    std::cout << std::endl;

    size_t out_i = 0;
    size_t in_i = 0;
    size_t total_frames_processed = 0;
    float* in_pos = signal->data();
    float* out_pos = out_signal->data();

    std::ofstream ftime("/tmp/ftim.csv");
    std::ofstream inp_file("/tmp/inp_file.csv");
    std::ofstream out_file("/tmp/out_file.csv");

    double overall_time = 0;
    while (in_i < in_signal_duration && out_i < out_signal_duration) {
        spx_uint32_t in_len = std::min(frame_size, in_signal_duration - in_i);
        spx_uint32_t out_len =
            std::min(static_cast<spx_uint32_t>(frame_size), static_cast<spx_uint32_t>(out_signal_duration - out_i));
        auto start_time = std::chrono::high_resolution_clock::now();
        err = speex_resampler_process_interleaved_float(resampler, in_pos, &in_len, out_pos, &out_len);
        const double duration =
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start_time)
                .count();
        const size_t navailable = out_len * channels;
        if (navailable > frame_size / 2) {
            const double rate = duration * 1e-9 / navailable * out_fs;
            ftime << rate << std::endl;
            time.push_back(rate);
            overall_time += duration;
        }

        for (size_t i = 0; i < in_len * channels; i++) {
            inp_file << in_pos[i] << std::endl;
        }
        for (size_t i = 0; i < out_len * channels; i++) {
            out_file << out_pos[i] << std::endl;
        }

        if (err != RESAMPLER_ERR_SUCCESS) {
            std::cerr << "Resampling error: " << err << std::endl;
            break;
        }

        out_i += out_len * channels;
        ;
        in_pos += in_len * channels;
        out_pos += out_len * channels;
        in_i += in_len * channels;
        total_frames_processed++;
    }
    double avg = 0;
    for (auto x : time) {
        avg += x;
    }
    avg = avg / time.size();
    double std = 0;
    for (auto x : time) {
        std += (x - avg) * (x - avg);
    }
    std = sqrt(std / time.size());

    // Calculate performance metrics
    std::cout << std::endl;
    std::cout << "Benchmark Results:" << std::endl;
    std::cout << "==================" << std::endl;
    std::cout << "Average ms in second: " << avg * 1000. << std::endl;
    std::cout << "STD ms in second: " << std * 1000. << std::endl;
    std::cout << "Overall time: " << overall_time * 1e-9 << std::endl;
    std::cout << "Realtimeness: " << overall_time * 1e-9 * static_cast<double>(out_fs) / static_cast<double>(out_i)
              << std::endl;

    // Clean up
    speex_resampler_destroy(resampler);

    ftime.close();
    inp_file.close();
    out_file.close();

    return 0;
}

// std::vector<double>&& do_speex(size_t in_fs, size_t out_fs, float* signal,
//     size_t in_signal_sz, float* out_signal, size_t *out_signal_sz, size_t frame_size)
// {
//     std::vector<double> time(*out_signal_sz);
//
//     src_t * src = src_open(SRC_PROFILE_MEDIUM, MONO, in_fs, out_fs);
//     size_t out_i = 0;
//     for (size_t i = 0; i < in_signal_sz && out_i < *out_signal_sz; i += frame_size) {
//         assert(src_push_samples(src, &signal[i], frame_size));
//         size_t navailable = 0;
//         do
//         {
//             auto start_time = std::chrono::high_resolution_clock::now();
//             err = speex_resampler_process_interleaved_float(resampler, in_pos, &in_len,
//                                               out_pos, &out_len);
//             const auto leftover = std::min(frame_size, out_signal->size() - out_i);
//             navailable = src_pop_samples(src, &(*out_signal)[out_i], leftover);
//             const double duration =
//             std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() -
//             start_time).count(); const size_t navailable = out_len * channels; if (navailable > frame_size / 2) {
//                 const double rate  = duration * 1e-9 / navailable * out_fs;
//                 ftime << rate << std::endl;
//                 time.push_back(rate);
//                 overall_time += duration;
//             }
//             out_i += navailable;
//         } while (navailable > 0 && out_i < out_signal_duration);
//     }
//
//     *out_signal_sz = out_i;
//     src_close(src);
// }
