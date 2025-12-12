SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})



add_executable(qr_benchmark
    qr.f90
    aux_print_matrix.f90
    qr_stats.f90
)

TARGET_LINK_LIBRARIES(qr_benchmark
        PRIVATE
            LPF::LPF
            LPF::BLAS
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )
