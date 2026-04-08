//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of LPF-BLAS, a BLAS Implementation for Half-Precision
   Copyright (C) 2025 Martin Koehler

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
   */

#include <nanobench.h>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
#include <iomanip>
#include "lpf_internal.h"
#include "fp16_helper.h"

// C-binding function declarations
extern "C" {
    lpf_ffloat16_t lpf_blas_hnrm2_fortran(lpf_blas_int_t *n,
                                           lpf_ffloat16_t *x,
                                           lpf_blas_int_t *incx);
    lpf_ffloat16_t lpf_blas_hnrm2_fp32_fortran(lpf_blas_int_t *n,
                                                 lpf_ffloat16_t *x,
                                                 lpf_blas_int_t *incx);
}

// Helper to convert fp16 to float
inline float fp16_to_float(lpf_ffloat16_t value) {
    union { _Float16 f16; uint16_t u16; } u;
    u.u16 = value.value;
    return static_cast<float>(u.f16);
}

// Helper to generate random fp16 data with controlled range to avoid overflow
void generate_random_fp16(lpf_ffloat16_t* data, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        // Generate random fp16 values in a reasonable range [-10, 10]
        float val = static_cast<float>(rand() % 2000 - 1000) / 100.0f;  // [-10, 10]
        union { float f32; uint32_t u32; } u;
        u.f32 = val;
        // Convert to fp16 by truncating
        uint16_t u16 = static_cast<uint16_t>(u.u32 >> 16);
        data[i].value = u16;
    }
}

// Reference implementation (compute in double to avoid overflow)
double reference_nrm2(const float* x, int n, int incx) {
    double ssq = 0.0;
    for (int i = 0; i < n; ++i) {
        double val = static_cast<double>(x[i * incx]);
        ssq += val * val;
    }
    return std::sqrt(ssq);
}

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    ankerl::nanobench::Bench b;
    b.performanceCounters(true);
    b.clockResolutionMultiple(10000000);
    b.minEpochTime(std::chrono::milliseconds(100));

    // Vector sizes to benchmark
    const std::vector<int> sizes = {32, 64, 128, 256, 512, 1024, 4096, 16384, 65536, 1048576};
    const std::vector<int> incs = {1, 2};

    // Fixed seed for reproducibility
    srand(12345);

    std::cout << "# hnrm2 / hnrm2_fp32 Benchmark" << std::endl;
    std::cout << "# nanobench version: " << ANKERL_NANOBENCH_VERSION_MAJOR << "."
              << ANKERL_NANOBENCH_VERSION_MINOR << "." << ANKERL_NANOBENCH_VERSION_PATCH << std::endl;
    std::cout << "# Performance counters (instructions, cpucycles, pagefaults)" << std::endl;
    std::cout << std::endl;

    for (int incx : incs) {
        for (int n : sizes) {
            // Allocate memory for the full array
            size_t array_size = static_cast<size_t>(n) * static_cast<size_t>(incx);
            std::vector<lpf_ffloat16_t> data(array_size);

            // Generate random data
            generate_random_fp16(data.data(), array_size);

            // Create view array with stride
            std::vector<lpf_ffloat16_t> x(n);
            for (int i = 0; i < n; ++i) {
                x[i] = data[i * incx];
            }

            // Compute reference result in double
            std::vector<float> x_float(n);
            for (int i = 0; i < n; ++i) {
                x_float[i] = fp16_to_float(x[i]);
            }
            double ref_result = reference_nrm2(x_float.data(), n, 1);

            // Clear previous results and benchmark hnrm2
            ankerl::nanobench::Bench hnrm2_bench;
            hnrm2_bench.performanceCounters(true);
            hnrm2_bench.clockResolutionMultiple(10000000);
            hnrm2_bench.minEpochTime(std::chrono::milliseconds(100));
            hnrm2_bench.run("hnrm2", [&x, n, incx]() {
                lpf_blas_int_t nn = n;
                lpf_blas_int_t incxx = incx;
                lpf_blas_hnrm2_fortran(&nn, x.data(), &incxx);
            });

            // Clear previous results and benchmark hnrm2_fp32
            ankerl::nanobench::Bench hnrm2_fp32_bench;
            hnrm2_fp32_bench.performanceCounters(true);
            hnrm2_fp32_bench.clockResolutionMultiple(10000000);
            hnrm2_fp32_bench.minEpochTime(std::chrono::milliseconds(100));
            hnrm2_fp32_bench.run("hnrm2_fp32", [&x, n, incx]() {
                lpf_blas_int_t nn = n;
                lpf_blas_int_t incxx = incx;
                lpf_blas_hnrm2_fp32_fortran(&nn, x.data(), &incxx);
            });

            // Get results
            auto results_hnrm2 = hnrm2_bench.results();
            auto results_hnrm2_fp32 = hnrm2_fp32_bench.results();

            if (!results_hnrm2.empty() && !results_hnrm2_fp32.empty()) {
                auto const& result_hnrm2 = results_hnrm2[0];
                auto const& result_hnrm2_fp32 = results_hnrm2_fp32[0];

                // Get median times
                double time_hnrm2 = result_hnrm2.median(ankerl::nanobench::Result::Measure::elapsed);
                double time_hnrm2_fp32 = result_hnrm2_fp32.median(ankerl::nanobench::Result::Measure::elapsed);

                // Output results
                std::cout << "n=" << std::setw(7) << n << ", incx=" << incx << ":" << std::endl;

                // hnrm2 result
                std::cout << "  hnrm2:           " << std::fixed << std::setprecision(2)
                          << (time_hnrm2 * 1e9) << " ns/op" << std::endl;

                // Print performance counters if available
                using Measure = ankerl::nanobench::Result::Measure;
                if (result_hnrm2.has(Measure::instructions)) {
                    double instructions = result_hnrm2.median(Measure::instructions);
                    double cpucycles = result_hnrm2.median(Measure::cpucycles);
                    std::cout << "    instructions:  " << static_cast<int64_t>(instructions) << std::endl;
                    std::cout << "    cpucycles:     " << static_cast<int64_t>(cpucycles) << std::endl;
                }

                // hnrm2_fp32 result
                std::cout << "  hnrm2_fp32:      " << std::fixed << std::setprecision(2)
                          << (time_hnrm2_fp32 * 1e9) << " ns/op" << std::endl;

                if (result_hnrm2_fp32.has(Measure::instructions)) {
                    double instructions = result_hnrm2_fp32.median(Measure::instructions);
                    double cpucycles = result_hnrm2_fp32.median(Measure::cpucycles);
                    std::cout << "    instructions:  " << static_cast<int64_t>(instructions) << std::endl;
                    std::cout << "    cpucycles:     " << static_cast<int64_t>(cpucycles) << std::endl;
                }

                std::cout << "  Reference:       " << std::fixed << std::setprecision(6)
                          << ref_result << std::endl;
                std::cout << std::endl;
            }
        }
    }

    return 0;
}
