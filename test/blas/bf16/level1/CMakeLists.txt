SET (SRC
    check_basum.f90
    check_baxpy.f90
    check_baxpby.f90
    check_bcopy.f90
    check_bdot.f90
    check_bsbdot.f90
    check_bnrm2.f90
    check_bnrm2_fp32.f90
    check_brot.f90
    check_brotg.f90
    check_brotm.f90
    check_brotmg.f90
    check_bscal.f90
    check_bswap.f90
    check_ibamax.f90
)

FOREACH( S IN LISTS SRC)
    cmake_path(GET S STEM LAST_ONLY check_name)
    ADD_EXECUTABLE(${check_name} ${S})

    TARGET_LINK_LIBRARIES(${check_name}
        PRIVATE
            LPF::lpf_blas
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )

    TARGET_INCLUDE_DIRECTORIES(${check_name}
        PRIVATE
            ${PROJECT_BINARY_DIR}/src
            ${PROJECT_SOURCE_DIR}/src
            $<$<BOOL:LPF_HAVE_BF16>:${PROJECT_SOURCE_DIR}/src/bf16>
    )

    ADD_TEST(NAME ${check_name} COMMAND ${check_name})

ENDFOREACH()
