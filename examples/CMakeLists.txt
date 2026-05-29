SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


SET (FP16_SRC
    fp16_demo_stats.f90
    fp16_demo1.f90
    fp16_demo_fraction.f90
)

SET (BF16_SRC
    bf16_demo_stats.f90
    bf16_demo1.f90
    bf16_demo_fraction.f90
)

SET(FP8_BLAS_SRC
    fp8_blas.f90
)

SET(FP8_BLAS_SRC
blas_cg.f90
)

SET(SRC)
SET(BLAS_SRC)

IF (LPF_HAVE_FP16)
    SET(SRC ${SRC} ${FP16_SRC})
    SET(BLAS_SRC ${BLAS_SRC} ${FP16_BLAS_SRC})
ENDIF()

IF (LPF_HAVE_BF16)
    SET(SRC ${SRC} ${BF16_SRC})
ENDIF()

IF (LPF_FP8_E4M3 AND LPF_FP8_E5M2)
    SET (BLAS_SRC ${BLAS_SRC} ${FP8_BLAS_SRC})
ENDIF()



FOREACH( S IN LISTS SRC)
    cmake_path(GET S STEM LAST_ONLY check_name)
    ADD_EXECUTABLE(${check_name} ${S})

    TARGET_LINK_LIBRARIES(${check_name}
        PRIVATE
            lpf
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )

ENDFOREACH()

if (LPF_BLAS)
    FOREACH( S IN LISTS BLAS_SRC)
    cmake_path(GET S STEM LAST_ONLY check_name)
    ADD_EXECUTABLE(${check_name} ${S})

    TARGET_LINK_LIBRARIES(${check_name}
        PRIVATE
            lpf
            lpf_blas
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )

ENDFOREACH()


endif()

if (LPF_EXAMPLES_TLAPACK)
    add_subdirectory(tlapack)
endif()
