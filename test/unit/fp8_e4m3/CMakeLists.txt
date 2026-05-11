# Shared test utilities module
ADD_LIBRARY(fp8_e4m3_test_utils STATIC lpf_fp8_e4m3_test_utils.f90)
TARGET_INCLUDE_DIRECTORIES(fp8_e4m3_test_utils
    PRIVATE
        ${PROJECT_BINARY_DIR}/src
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}/src/fortran-modules
)
TARGET_LINK_LIBRARIES(fp8_e4m3_test_utils
    PRIVATE
        lpf
        $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
)

# List of unit test executables
SET(FP8_E4M3_UNIT_TESTS
    test_fp8_e4m3_constructors
    test_fp8_e4m3_assignment
    test_fp8_e4m3_arithmetic
    test_fp8_e4m3_comparison
    test_fp8_e4m3_math_one_arg
    test_fp8_e4m3_math_two_arg
    test_fp8_e4m3_trig_deg
    test_fp8_e4m3_intrinsics
    test_fp8_e4m3_minmax
    test_fp8_e4m3_reductions
    test_fp8_e4m3_io
)

FOREACH(test_name IN LISTS FP8_E4M3_UNIT_TESTS)
    ADD_EXECUTABLE(${test_name} ${test_name}.f90)
    TARGET_LINK_LIBRARIES(${test_name}
        PRIVATE
            lpf
            fp8_e4m3_test_utils
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )
    TARGET_INCLUDE_DIRECTORIES(${test_name}
        PRIVATE
            ${PROJECT_BINARY_DIR}/src
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_BINARY_DIR}/src/fortran-modules
    )
    ADD_TEST(NAME UNIT-FP8_E4M3-${test_name} COMMAND ${test_name})
ENDFOREACH()