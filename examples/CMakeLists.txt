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

SET(SRC)
IF (LPF_HAVE_FP16)
    SET(SRC ${SRC} ${FP16_SRC})
ENDIF()

IF (LPF_HAVE_BF16)
    SET(SRC ${SRC} ${BF16_SRC})
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
    add_subdirectory(qr)
endif()

# ADD_SUBDIRECTORY(benchmark)
