SET(CMAKE_Fortran_STANDARD 08)
SET(SYSINFO_SRC sysinfo/${CMAKE_SYSTEM_PROCESSOR}.c)


SET(FP16_SRC
    fp16/fp16.F90
    fp16/fp16_math_generic.f90
    fp16/fp16_math_one_arg.f90
    fp16/fp16_math_two_args.f90
    fp16/fp16_math_trig.f90
    fp16/fp16_maxval.f90
    fp16/fp16_maxloc.f90
    fp16/fp16_minloc.f90
    fp16/fp16_minval.f90
    fp16/fp16_helper.c
    fp16/fp16_math_one_arg.c
    fp16/fp16_math_two_args.c
    fp16/fp16_math_generic.c
    fp16/fp16_math_trig.c
)

SET (BF16_SRC
    bf16/bf16.F90
    bf16/bf16_math_generic.f90
    bf16/bf16_math_one_arg.f90
    bf16/bf16_math_two_args.f90
    bf16/bf16_math_trig.f90
    bf16/bf16_maxval.f90
    bf16/bf16_maxloc.f90
    bf16/bf16_minloc.f90
    bf16/bf16_minval.f90
    bf16/bf16_helper.c
    bf16/bf16_math_one_arg.c
    bf16/bf16_math_two_args.c
    bf16/bf16_math_generic.c
    bf16/bf16_math_trig.c
)

SET(FP8_E4M3_SRC
    fp8_e4m3/emu/classify.c
    fp8_e4m3/emu/conversion.c
    fp8_e4m3/emu/base_arithmetic.c
    fp8_e4m3/emu/debug.c
    fp8_e4m3/emu/fp8_e4m3_math_lut.c
    fp8_e4m3/emu/fp8_e4m3_abs.c
    fp8_e4m3/fp8_e4m3.F90
    fp8_e4m3/fp8_e4m3_generic.c
    fp8_e4m3/fp8_e4m3_helper.c
    fp8_e4m3/fp8_e4m3_math_generic.f90
    fp8_e4m3/fp8_e4m3_math_one_arg.c
    fp8_e4m3/fp8_e4m3_math_one_arg.f90
    fp8_e4m3/fp8_e4m3_math_two_args.c
    fp8_e4m3/fp8_e4m3_math_two_args.f90
    fp8_e4m3/fp8_e4m3_math_trig.c
    fp8_e4m3/fp8_e4m3_math_trig.f90
    fp8_e4m3/fp8_e4m3_maxval.f90
    fp8_e4m3/fp8_e4m3_maxloc.f90
    fp8_e4m3/fp8_e4m3_minloc.f90
    fp8_e4m3/fp8_e4m3_minval.f90
)

SET(FP8_E5M2_SRC
    fp8_e5m2/emu/classify.c
    fp8_e5m2/emu/conversion.c
    fp8_e5m2/emu/base_arithmetic.c
    fp8_e5m2/emu/debug.c
    fp8_e5m2/emu/abs.c
    fp8_e5m2/emu/fp8_e5m2_math_lut.c
    fp8_e5m2/fp8_e5m2.F90
    fp8_e5m2/fp8_e5m2_generic.c
    fp8_e5m2/fp8_e5m2_helper.c
    fp8_e5m2/fp8_e5m2_math_generic.f90
    fp8_e5m2/fp8_e5m2_math_one_arg.c
    fp8_e5m2/fp8_e5m2_math_one_arg.f90
    fp8_e5m2/fp8_e5m2_math_two_args.c
    fp8_e5m2/fp8_e5m2_math_two_args.f90
    fp8_e5m2/fp8_e5m2_math_trig.c
    fp8_e5m2/fp8_e5m2_math_trig.f90
    fp8_e5m2/fp8_e5m2_maxval.f90
    fp8_e5m2/fp8_e5m2_maxloc.f90
    fp8_e5m2/fp8_e5m2_minloc.f90
    fp8_e5m2/fp8_e5m2_minval.f90
)

SET(PUBLIC_HEADERS
    $<$<BOOL:${LPF_HAVE_FP16}>:${CMAKE_CURRENT_SOURCE_DIR}/cxx/lpf_fp16.hpp>
    $<$<BOOL:${LPF_HAVE_FP16}>:${CMAKE_CURRENT_SOURCE_DIR}/cxx/lpf_fp16_tlapack.hpp>
    $<$<BOOL:${LPF_HAVE_BF16}>:${CMAKE_CURRENT_SOURCE_DIR}/cxx/lpf_bf16.hpp>
    $<$<BOOL:${LPF_HAVE_BF16}>:${CMAKE_CURRENT_SOURCE_DIR}/cxx/lpf_bf16_tlapack.hpp>
    $<$<BOOL:${LPF_FP8_E4M3}>:${CMAKE_CURRENT_SOURCE_DIR}/fp8_e4m3/emu/lpf_fp8_e4m3_emu.h>
    $<$<BOOL:${LPF_FP8_E5M2}>:${CMAKE_CURRENT_SOURCE_DIR}/fp8_e5m2/emu/lpf_fp8_e5m2_emu.h>
)


ADD_LIBRARY(lpf
    $<$<BOOL:${LPF_HAVE_FP16}>:${FP16_SRC}>
    $<$<BOOL:${LPF_HAVE_BF16}>:${BF16_SRC}>
    $<$<BOOL:${LPF_FP8_E4M3}>:${FP8_E4M3_SRC}>
    $<$<BOOL:${LPF_FP8_E5M2}>:${FP8_E5M2_SRC}>
    ${PUBLIC_HEADERS}
    aux/get_wtime.c
    aux/xerbla.c
    aux/xerbla.f90
    init.c
    types.F90
    ${SYSINFO_SRC}
)

TARGET_LINK_LIBRARIES(lpf
    PUBLIC
        $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
)


TARGET_INCLUDE_DIRECTORIES(lpf
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cxx>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fp8_e4m3/emu>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fp8_e5m2/emu>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules>
        $<INSTALL_INTERFACE:include/lpf>
)

TARGET_COMPILE_DEFINITIONS(lpf PUBLIC
    $<$<BOOL:${LPF_HAVE_FP16}>:-DLPF_HAVE_FP16>
    $<$<BOOL:${LPF_HAVE_BF16}>:-DLPF_HAVE_BF16>
)

if(LPF_USE_IPO)
  set_property(TARGET lpf PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

SET_TARGET_PROPERTIES(lpf PROPERTIES PUBLIC_HEADER
    "${PUBLIC_HEADERS}"
)
SET_PROPERTY(TARGET lpf PROPERTY VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
SET_PROPERTY(TARGET lpf PROPERTY SOVERSION ${PROJECT_VERSION_MAJOR})
SET_TARGET_PROPERTIES(lpf PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules/" DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR} )

ADD_LIBRARY(LPF::lpf ALIAS lpf)

INSTALL(TARGETS lpf
        EXPORT  lpf_targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lpf
        )

if(LPF_BLAS)
    ADD_SUBDIRECTORY(blas)
endif()

enable_testing()
#
# Tests
#
if (LPF_FP8_E4M3)
    add_executable(test_fp8_e4m3 fp8_e4m3/emu/tests.c)
    target_link_libraries(test_fp8_e4m3 PRIVATE lpf)
    add_test(NAME test_fp8_e4m3 COMMAND test_fp8_e4m3 )

    # Lookup table generator (developer tool)
    add_executable(gen_fp8_e4m3_lut fp8_e4m3/emu/gen_lookup_tables.c
        fp8_e4m3/emu/conversion.c fp8_e4m3/emu/classify.c)
    target_include_directories(gen_fp8_e4m3_lut PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(gen_fp8_e4m3_lut PRIVATE m)
endif()

if (LPF_FP8_E5M2)
    add_executable(test_fp8_e5m2 fp8_e5m2/emu/tests.c)
    target_link_libraries(test_fp8_e5m2 PRIVATE lpf)
    add_test(NAME test_fp8_e5m2 COMMAND test_fp8_e5m2 )

    # Lookup table generator (developer tool)
    add_executable(gen_fp8_e5m2_lut fp8_e5m2/emu/gen_lookup_tables.c
        fp8_e5m2/emu/conversion.c fp8_e5m2/emu/classify.c)
    target_include_directories(gen_fp8_e5m2_lut PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(gen_fp8_e5m2_lut PRIVATE m)
endif()
