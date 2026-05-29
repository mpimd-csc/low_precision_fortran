# Detect Externals
IF (LPF_INTEGER8)
    SET(BLA_SIZEOF_INTEGER 8)
ELSE()
    SET(BLA_SIZEOF_INTEGER 4)
ENDIF()
FIND_PACKAGE(BLAS)
FIND_PACKAGE(LAPACK)

# Search for special in routines in BLAS
IF (BLAS_FOUND AND LAPACK_FOUND)
SET(_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
SET(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
check_function_exists(cblas_hgemm LPF_BLAS_USE_HGEMM)
check_function_exists(cblas_gemm_f16f16f32 LPF_BLAS_USE_GEMM_F16F16F32)
check_function_exists(cblas_gemm_bf16bf16f32 LPF_BLAS_USE_GEMM_BF16BF16F32)
check_function_exists(aocl_gemm_bf16bf16f32of32 LPF_BLAS_USE_AMD_GEMM_BF16BF16F32)
check_function_exists(mkl_get_version LPF_BLAS_IS_MKL)
check_function_exists(AOCL_BLIS_set_progress LPF_BLAS_IS_AOCL)
SET(CMAKE_REQUIRED_LIBRARIES ${_CMAKE_REQUIRED_LIBRARIES})

IF (LPF_BLAS_IS_MKL)
    FIND_PATH(MKL_INCLUDE_DIR mkl.h
        HINTS
        "$ENV{MKLROOT}/include"
        "${MKLROOT}/include"
        "$ENV{MKL_ROOT}/include"
        "${MKL_ROOT}/include"
    )
    if (NOT MKL_INCLUDE_DIR)
        message(FATAL_ERROR "Inte (one)MKL detected but mkl.h not found.")
    endif()
    message(STATUS "MKL include directory: ${MKL_INCLUDE_DIR}")
    set_property(TARGET BLAS::BLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
ENDIF()
ENDIF()

SET (SRC )
if (LPF_HAVE_FP16)
    ADD_SUBDIRECTORY(fp16/level1)
    ADD_SUBDIRECTORY(fp16/level2)
    ADD_SUBDIRECTORY(fp16/level3)
    ADD_SUBDIRECTORY(fp16/aux)
    SET(SRC ${SRC}
        fp16/scaling/scale_diag_twoside.f90
        fp16/scaling/scale_diag_right.f90
        fp16/lpf_blas_fp16.F90

        fp16/lapack/cholesky/potf2.f90
        fp16/lapack/cholesky/potrf.f90
        fp16/lapack/cholesky/potrf2.f90
        fp16/lapack/cholesky/potrfp.f90

       fp16/lapack/lange.f90
       fp16/lpf_lapack_fp16.F90


    )
endif()

if (LPF_HAVE_BF16)
    ADD_SUBDIRECTORY(bf16/level1)
    ADD_SUBDIRECTORY(bf16/level2)
    ADD_SUBDIRECTORY(bf16/level3)
    ADD_SUBDIRECTORY(bf16/aux)
    SET(SRC ${SRC}
        bf16/scaling/scale_diag_twoside.f90
        bf16/scaling/scale_diag_right.f90
        bf16/lpf_blas_bf16.F90

        bf16/lapack/cholesky/potf2.f90
        bf16/lapack/cholesky/potrf.f90
        bf16/lapack/cholesky/potrf2.f90
        bf16/lapack/cholesky/potrfp.f90

        bf16/lapack/lange.f90
        bf16/lpf_lapack_bf16.F90
    )
endif()


SET(SRC ${SRC}
    highprecision/scale_diag_twoside.f90
    highprecision/scale_diag_right.f90
)

SET(SRC_GENERIC
    generic/lpf_blas_generic.F90
    generic/level1/asum.F90
    generic/level1/axpy.F90
    generic/level1/axpby.F90
    generic/level1/copy.F90
    generic/level1/dot.F90
    generic/level1/nrm2.F90
    generic/level1/nrm2_fp32.F90
    generic/level1/rot.F90
    generic/level1/rotg.F90
    generic/level1/rotmg.F90
    generic/level1/rotm.F90
    generic/level1/iamax.F90
    generic/level1/sbdot.F90
    generic/level1/scal.F90
    generic/level1/swap.F90

    generic/level2/gbmv.F90
    generic/level2/gemv.F90
    generic/level2/ger.F90
    generic/level2/sbmv.F90
    generic/level2/spmv.F90
    generic/level2/spr.F90
    generic/level2/spr2.F90
    generic/level2/symv.F90
    generic/level2/tbmv.F90
    generic/level2/tbsv.F90
    generic/level2/tpmv.F90
    generic/level2/tpsv.F90
    generic/level2/trmv.F90
    generic/level2/trsv.F90

    generic/level3/gemm.F90
    generic/level3/gemm_fp32.F90
    generic/level3/gemmtr.F90
    generic/level3/symm.F90
    generic/level3/syr2.F90
    generic/level3/syr2k.F90
    generic/level3/syrk.F90
    generic/level3/trmm.F90
    generic/level3/trsm.F90

)

if (LPF_FP8_E5M2)
    ADD_LIBRARY(lpf_blas_fp8_e5m2 OBJECT ${SRC_GENERIC})
    target_compile_definitions(lpf_blas_fp8_e5m2 PRIVATE LPF_FP8_E5M2 DT=fp8_e5m2 )
    set_target_properties(lpf_blas_fp8_e5m2 PROPERTIES  Fortran_PREPROCESS  ON )
    target_link_libraries(lpf_blas_fp8_e5m2 PRIVATE lpf)
    set_target_properties(lpf_blas_fp8_e5m2 PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
endif()

if (LPF_FP8_E4M3)
    ADD_LIBRARY(lpf_blas_fp8_e4m3 OBJECT ${SRC_GENERIC})
    target_compile_definitions(lpf_blas_fp8_e4m3 PRIVATE LPF_FP8_E4M3 DT=fp8_e4m3 )
    set_target_properties(lpf_blas_fp8_e4m3 PROPERTIES  Fortran_PREPROCESS  ON )
    target_link_libraries(lpf_blas_fp8_e4m3 PRIVATE lpf)
    set_target_properties(lpf_blas_fp8_e4m3 PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
endif()




ADD_LIBRARY(lpf_blas
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level1>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level2>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level3>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_aux>>

    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level1>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level2>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level3>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_aux>>

    $<$<BOOL:${LPF_FP8_E5M2}>:$<TARGET_OBJECTS:lpf_blas_fp8_e5m2>>
    $<$<BOOL:${LPF_FP8_E4M3}>:$<TARGET_OBJECTS:lpf_blas_fp8_e4m3>>

    ${SRC}
)

TARGET_LINK_LIBRARIES(lpf_blas
    PRIVATE
    $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    PUBLIC
    LPF::lpf
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)

TARGET_INCLUDE_DIRECTORIES(lpf_blas
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules>
    $<INSTALL_INTERFACE:include/lpf>
)

TARGET_COMPILE_DEFINITIONS(lpf_blas PRIVATE $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>)

if(LPF_USE_IPO)
    set_property(TARGET lpf_blas PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

SET_PROPERTY(TARGET lpf_blas PROPERTY VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
SET_PROPERTY(TARGET lpf_blas PROPERTY SOVERSION ${PROJECT_VERSION_MAJOR})
SET_TARGET_PROPERTIES(lpf_blas PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules/" DESTINATION ${CMAKE_INSTALL_FULL_MODULEDIR} )

ADD_LIBRARY(LPF::lpf_blas ALIAS lpf_blas)

INSTALL(TARGETS lpf_blas
    EXPORT  lpf_targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)



