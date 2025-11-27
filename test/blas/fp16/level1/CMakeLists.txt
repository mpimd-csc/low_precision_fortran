SET (SRC
    check_hasum.c
    check_haxpy.c
    check_haxpby.c
    check_hcopy.c
    check_hdot.c
    check_hshdot.c
    check_hnrm2.c
    check_hnrm2_fp32.c
    check_hrot.c
    check_hrotg.c
    check_hrotm.c
    check_hrotmg.c
    check_hscal.c
    check_hswap.c
    check_ihamax.c
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
            $<$<BOOL:LPF_HAVE_FP16>:${PROJECT_SOURCE_DIR}/src/fp16>
    )

    ADD_TEST(NAME ${check_name} COMMAND ${check_name})

ENDFOREACH()

# ADD_EXECUTABLE(fp16_level1
#     main.c
#     header.c
#     htest.c
#     htest1.c
#     itest1.c
#     check0.c
#     check1.c
#     check_hasum.c
#     check_hnrm2.c
#     check_hscal.c
#     check_ihamax.c
# )
#
# TARGET_LINK_LIBRARIES(fp16_level1
#     PRIVATE
#         hpblas
#         $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
# )
#
# TARGET_INCLUDE_DIRECTORIES(fp16_level1
#     PRIVATE
#         ${PROJECT_BINARY_DIR}/src
#         ${PROJECT_SOURCE_DIR}/src
# )
#
# ADD_TEST(NAME FP16_LEVEL1 COMMAND fp16_level1)
#
