# Shared test utilities module
ADD_LIBRARY(fp16_test_utils STATIC lpf_fp16_test_utils.f90)
add_coverage(fp16_test_utils)
TARGET_INCLUDE_DIRECTORIES(fp16_test_utils
    PRIVATE
        ${PROJECT_BINARY_DIR}/src
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}/src/fortran-modules
)
TARGET_LINK_LIBRARIES(fp16_test_utils
    PRIVATE
        lpf
        $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
)

# List of unit test executables
SET(FP16_UNIT_TESTS
    test_fp16_constructors
    test_fp16_assignment
    test_fp16_arithmetic
    test_fp16_comparison
    test_fp16_math_one_arg
    test_fp16_math_two_arg
    test_fp16_trig_deg
    test_fp16_intrinsics
    test_fp16_minmax
    test_fp16_reductions
    test_fp16_io
    test_fp16_matmul
    test_fp16_dot_product
)

FOREACH(test_name IN LISTS FP16_UNIT_TESTS)
    ADD_EXECUTABLE(${test_name} ${test_name}.f90)
    add_coverage(${test_name})
    TARGET_LINK_LIBRARIES(${test_name}
        PRIVATE
            lpf
            fp16_test_utils
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )
    TARGET_INCLUDE_DIRECTORIES(${test_name}
        PRIVATE
            ${PROJECT_BINARY_DIR}/src
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_BINARY_DIR}/src/fortran-modules
    )
    ADD_TEST(NAME UNIT-FP16-${test_name} COMMAND ${test_name})
ENDFOREACH()
