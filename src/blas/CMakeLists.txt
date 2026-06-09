FIND_PACKAGE(BLAS)
FIND_PACKAGE(LAPACK)

# Search for special in routines in BLAS
IF (BLAS_FOUND AND LAPACK_FOUND)
    SET(_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
SET(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
include(CheckCSourceRuns)

# Store the current linker status for the test
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})

# Test source code: Passes 64-bit integer via pointer to BLAS
check_c_source_runs("
#include <stdint.h>
#include <stdio.h>
double ddot_(const int64_t* n, const double* dx, const int64_t* incx,  const double* dy, const int64_t* incy  );

int main() {
    int32_t nx[3] = { 1, 1, 1 };
    double x[1] = {1.0};
    double y[1] = {1.0};
    int64_t incx = 0;
    int64_t incy = 0;

    double d = ddot_((int64_t *)&nx[1], x, &incx, y, &incy);
    printf(\"d = %lg\", d);
    if ( d <= 1.0) {
        return 1;
    }
    return 0;
}
" BLAS_IS_ILP64)

if(BLAS_IS_ILP64)
    message(STATUS "BLAS Integer Model: ILP64.")
else()
    message(STATUS "BLAS Integer Model: LP64.")
endif()

check_function_exists(cblas_hgemm LPF_BLAS_USE_HGEMM)
check_function_exists(cblas_gemm_f16f16f32 LPF_BLAS_USE_GEMM_F16F16F32)
check_function_exists(cblas_gemm_bf16bf16f32 LPF_BLAS_USE_GEMM_BF16BF16F32)
check_function_exists(cblas_gemm_e5m2e5m2f32 LPF_BLAS_USE_GEMM_E5M2_F32)
check_function_exists(cblas_gemm_e4m3e4m3f32 LPF_BLAS_USE_GEMM_E4M3_F32)
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
    SET (SRC ${SRC}
        fp16/level1/hasum.c
        fp16/level1/haxpby.c
        fp16/level1/haxpy.c
        fp16/level1/hcopy.c
        fp16/level1/hdot.c
        fp16/level1/hnrm2.c
        fp16/level1/hnrm2_fp32.c
        fp16/level1/hrot.c
        fp16/level1/hrotg.c
        fp16/level1/hrotm.c
        fp16/level1/hrotmg.c
        fp16/level1/hscal.c
        fp16/level1/hshdot.c
        fp16/level1/hswap.c
        fp16/level1/ihamax.c

        fp16/level2/hgbmv.c
        fp16/level2/hgemv.c
        fp16/level2/hger.c
        fp16/level2/hsbmv.c
        fp16/level2/hspmv.c
        fp16/level2/hspr2.c
        fp16/level2/hspr.c
        fp16/level2/hsymv.c
        fp16/level2/hsyr2.c
        fp16/level2/hsyr.c
        fp16/level2/htbmv.c
        fp16/level2/htbsv.c
        fp16/level2/htpmv.c
        fp16/level2/htpsv.c
        fp16/level2/htrmv.c
        fp16/level2/htrsv.c

        fp16/level3/hgemm.c
        fp16/level3/hgemm_fp32.c
        fp16/level3/hgemmtr.c
        fp16/level3/hsymm.c
        fp16/level3/hsyr2k.c
        fp16/level3/hsyrk.c
        fp16/level3/htrmm.c
        fp16/level3/htrsm.c

        fp16/aux/h2scopy.c
        fp16/aux/h2slacpy.c
        fp16/aux/hlacpy.c
        fp16/aux/hlaset.c
        fp16/aux/s2hcopy.c
        fp16/aux/s2hlacpy.c
    )
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
    SET( SRC ${SRC}
        bf16/level1/basum.c
        bf16/level1/baxpby.c
        bf16/level1/baxpy.c
        bf16/level1/bcopy.c
        bf16/level1/bdot.c
        bf16/level1/bnrm2.c
        bf16/level1/bnrm2_fp32.c
        bf16/level1/brot.c
        bf16/level1/brotg.c
        bf16/level1/brotm.c
        bf16/level1/brotmg.c
        bf16/level1/bscal.c
        bf16/level1/bsbdot.c
        bf16/level1/bswap.c
        bf16/level1/ibamax.c

        bf16/level2/bgbmv.c
        bf16/level2/bgemv.c
        bf16/level2/bger.c
        bf16/level2/bsbmv.c
        bf16/level2/bspmv.c
        bf16/level2/bspr2.c
        bf16/level2/bspr.c
        bf16/level2/bsymv.c
        bf16/level2/bsyr2.c
        bf16/level2/bsyr.c
        bf16/level2/btbmv.c
        bf16/level2/btbsv.c
        bf16/level2/btpmv.c
        bf16/level2/btpsv.c
        bf16/level2/btrmv.c
        bf16/level2/btrsv.c

        bf16/level3/bgemm.c
        bf16/level3/bgemm_fp32.c
        bf16/level3/bgemmtr.c
        bf16/level3/bsymm.c
        bf16/level3/bsyr2k.c
        bf16/level3/bsyrk.c
        bf16/level3/btrmm.c
        bf16/level3/btrsm.c
        bf16/aux/b2scopy.c
        bf16/aux/b2slacpy.c
        bf16/aux/blacpy.c
        bf16/aux/blaset.c
        bf16/aux/s2bcopy.c
        bf16/aux/s2blacpy.c
    )

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

    ADD_LIBRARY(lpf_blas_fp8_e5m2
        OBJECT
        ${SRC_GENERIC}
        $<$<BOOL:${LPF_BLAS_USE_GEMM_E5M2_F32}>:generic/interfaces/cblas_gemm_e5m2_fp32.c>
    )
    target_compile_definitions(lpf_blas_fp8_e5m2 PRIVATE
        LPF_FP8_E5M2
        DT=fp8_e5m2
        $<$<BOOL:${BLAS_IS_ILP64}>:BLAS_IS_ILP64=1>
        $<$<BOOL:${LPF_BLAS_USE_GEMM_E5M2_F32}>:LPF_BLAS_USE_GEMM_E5M2_F32>
    )
    set_target_properties(lpf_blas_fp8_e5m2 PROPERTIES  Fortran_PREPROCESS  ON )
    target_link_libraries(lpf_blas_fp8_e5m2 PRIVATE lpf)
    set_target_properties(lpf_blas_fp8_e5m2 PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
    add_coverage(lpf_blas_fp8_e5m2)

endif()

if (LPF_FP8_E4M3)
    ADD_LIBRARY(lpf_blas_fp8_e4m3 OBJECT
        ${SRC_GENERIC}
        $<$<BOOL:${LPF_BLAS_USE_GEMM_E4M3_F32}>:generic/interfaces/cblas_gemm_e4m3_fp32.c>
    )
    target_compile_definitions(lpf_blas_fp8_e4m3 PRIVATE
        LPF_FP8_E4M3
        DT=fp8_e4m3
        $<$<BOOL:${BLAS_IS_ILP64}>:BLAS_IS_ILP64=1>
        $<$<BOOL:${LPF_BLAS_USE_GEMM_E4M3_F32}>:LPF_BLAS_USE_GEMM_E4M3_F32>
    )
    set_target_properties(lpf_blas_fp8_e4m3 PROPERTIES  Fortran_PREPROCESS  ON )
    target_link_libraries(lpf_blas_fp8_e4m3 PRIVATE lpf)
    set_target_properties(lpf_blas_fp8_e4m3 PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran-modules")
    add_coverage(lpf_blas_fp8_e4m3)
endif()




ADD_LIBRARY(lpf_blas
    $<$<BOOL:${LPF_FP8_E5M2}>:$<TARGET_OBJECTS:lpf_blas_fp8_e5m2>>
    $<$<BOOL:${LPF_FP8_E4M3}>:$<TARGET_OBJECTS:lpf_blas_fp8_e4m3>>

    ${SRC}
)


#
# Fix strange optimization behavior
#
IF ( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    MESSAGE(STATUS "Remove too strong optimization from nrm2_fp32 computation on x86_64.")

    set_source_files_properties(fp16/level1/hnrm2_fp32.c PROPERTIES COMPILE_OPTIONS "-mno-avx512fp16;-mno-avx512bf16")
    set_source_files_properties(bf16/level1/bnrm2_fp32.c PROPERTIES COMPILE_OPTIONS "-mno-avx512fp16;-mno-avx512bf16")
ENDIF()





add_coverage(lpf_blas)

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
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/fp16>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/bf16>
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran-modules>
    $<INSTALL_INTERFACE:include/lpf>
)

TARGET_COMPILE_DEFINITIONS(lpf_blas PRIVATE
    $<$<BOOL:${LPF_BLAS_IS_MKL}>:-DLPF_BLAS_IS_MKL>
    $<$<BOOL:${LPF_BLAS_USE_GEMM_BF16BF16F32}>:-DLPF_BLAS_USE_GEMM_BF16BF16F32>
    $<$<BOOL:${LPF_BLAS_USE_HGEMM}>:-DLPF_BLAS_USE_HGEMM>
    $<$<BOOL:${LPF_BLAS_USE_GEMM_F16F16F32}>:-DLPF_BLAS_USE_GEMM_F16F16F32>
    $<$<BOOL:${BLAS_IS_ILP64}>:BLAS_IS_ILP64=1>
)

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



