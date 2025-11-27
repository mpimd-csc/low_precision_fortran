SET (SRC
    check_basum.c
    check_baxpy.c
    check_baxpby.c
    check_bcopy.c
    check_bdot.c
    check_bsbdot.c
    check_bnrm2.c
    check_bnrm2_fp32.c
    check_brot.c
    check_brotg.c
    check_brotm.c
    check_brotmg.c
    check_bscal.c
    check_bswap.c
    check_ibamax.c
)

FOREACH( S IN LISTS SRC)
    cmake_path(GET S STEM LAST_ONLY check_name)
    ADD_EXECUTABLE(${check_name} ${S})

    TARGET_LINK_LIBRARIES(${check_name}
        PRIVATE
            LPF::BLAS
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
