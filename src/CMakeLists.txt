SET(CMAKE_Fortran_STANDARD 08)
SET(FP16_SRC
    fp16/fp16.f90
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
    bf16/bf16.f90
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

ADD_LIBRARY(lpf
    $<$<BOOL:${LPF_HAVE_FP16}>:${FP16_SRC}>
    $<$<BOOL:${LPF_HAVE_BF16}>:${BF16_SRC}>
)

TARGET_LINK_LIBRARIES(lpf
    PRIVATE
        $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
)

TARGET_INCLUDE_DIRECTORIES(lpf
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules>
        $<INSTALL_INTERFACE:include/lpf>
)

if(LPF_USE_IPO)
  set_property(TARGET lpf PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

SET_PROPERTY(TARGET lpf PROPERTY VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
SET_PROPERTY(TARGET lpf PROPERTY SOVERSION ${PROJECT_VERSION_MAJOR})
SET_TARGET_PROPERTIES(lpf PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules/" DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR} )

ADD_LIBRARY(LPF::LPF ALIAS lpf)

INSTALL(TARGETS lpf
        EXPORT  lpf_targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )

