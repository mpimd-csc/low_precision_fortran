SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

ENABLE_LANGUAGE(CXX)
set(CMAKE_CXX_STANDARD 17)

# Fetch nanobench
include(FetchContent)
FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG v4.3.11
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nanobench)

add_executable(qr_benchmark
    qr.f90
    aux_print_matrix.f90
    qr_stats.f90
    qr_calls.f90
)

add_executable(bench_hnrm2
    bench_hnrm2.cpp
)
target_include_directories(bench_hnrm2 PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/fp16
)
target_link_libraries(bench_hnrm2 PRIVATE LPF::LPF LPF::BLAS nanobench)

TARGET_LINK_LIBRARIES(qr_benchmark
        PRIVATE
            LPF::LPF
            LPF::BLAS
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )
