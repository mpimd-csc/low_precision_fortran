# Shared test utilities module
ADD_LIBRARY(bf16_test_utils STATIC lpf_bf16_test_utils.f90)
TARGET_INCLUDE_DIRECTORIES(bf16_test_utils
    PRIVATE
        ${PROJECT_BINARY_DIR}/src
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}/src/fortran-modules
)
TARGET_LINK_LIBRARIES(bf16_test_utils
    PRIVATE
        lpf
        $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
)

# List of unit test executables
SET(BF16_UNIT_TESTS
    test_bf16_constructors
    test_bf16_assignment
    test_bf16_arithmetic
    test_bf16_comparison
    test_bf16_math_one_arg
    test_bf16_math_two_arg
    test_bf16_trig_deg
    test_bf16_intrinsics
    test_bf16_minmax
    test_bf16_reductions
    test_bf16_io
)

FOREACH(test_name IN LISTS BF16_UNIT_TESTS)
    ADD_EXECUTABLE(${test_name} ${test_name}.f90)
    TARGET_LINK_LIBRARIES(${test_name}
        PRIVATE
            lpf
            bf16_test_utils
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )
    TARGET_INCLUDE_DIRECTORIES(${test_name}
        PRIVATE
            ${PROJECT_BINARY_DIR}/src
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_BINARY_DIR}/src/fortran-modules
    )
    ADD_TEST(NAME UNIT-BF16-${test_name} COMMAND ${test_name})
ENDFOREACH()