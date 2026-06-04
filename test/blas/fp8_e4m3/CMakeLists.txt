# Create the utility library
SET(SUFFIX "_fp8_e4m3")
ADD_LIBRARY(generic_test_utils${SUFFIX} STATIC ../generic_src/generic_test_utils.F90)
target_compile_definitions(generic_test_utils${SUFFIX} PRIVATE LPF_FP8_E4M3 DT=fp8_e4m3 )
TARGET_LINK_LIBRARIES(generic_test_utils${SUFFIX} PRIVATE lpf lpf_blas)

# List of test programs
SET(GENERIC_LEVEL1_TESTS
    test_level1_asum
    test_level1_axpby
    test_level1_axpy
    test_level1_copy
    test_level1_dot
    test_level1_iamax
    test_level1_nrm2
    test_level1_nrm2_fp32
    test_level1_rot
    test_level1_rotg
    test_level1_rotm
    test_level1_rotmg
    test_level1_scal
    test_level1_swap
)

SET(GENERIC_LEVEL2_TESTS
    test_level2_gbmv
    test_level2_gemv
    test_level2_ger
    test_level2_sbmv
    test_level2_spmv
    test_level2_spr
    test_level2_spr2
    test_level2_symv
    test_level2_tbmv
    test_level2_tbsv
    test_level2_tpmv
    test_level2_tpsv
    test_level2_trmv
    test_level2_trsv
)

SET(GENERIC_LEVEL3_TESTS
    test_level3_gemm
    test_level3_gemm_fp32
    test_level3_symm
    test_level3_syrk
    test_level3_syr2k
    test_level3_trmm
    test_level3_trsm
    test_level3_gemmtr
)

SET(GENERIC_TESTS ${GENERIC_LEVEL1_TESTS} ${GENERIC_LEVEL2_TESTS} ${GENERIC_LEVEL3_TESTS})

FOREACH(TEST_NAME ${GENERIC_TESTS})
    ADD_EXECUTABLE(${TEST_NAME}${SUFFIX} ../generic_src/${TEST_NAME}.F90)
    add_coverage(${TEST_NAME}${SUFFIX})
    TARGET_LINK_LIBRARIES(${TEST_NAME}${SUFFIX} PRIVATE lpf lpf_blas generic_test_utils${SUFFIX})
    target_compile_definitions(${TEST_NAME}${SUFFIX} PRIVATE LPF_FP8_E4M3 DT=fp8_e4m3 )
    set_target_properties(${TEST_NAME}${SUFFIX} PROPERTIES  Fortran_PREPROCESS  ON )
    ADD_TEST(NAME "UNIT-BLAS-${TEST_NAME}${SUFFIX}" COMMAND ${TEST_NAME}${SUFFIX})
ENDFOREACH()

