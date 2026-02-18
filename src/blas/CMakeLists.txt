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
        fp16/lapack/cholqr/gecholqr.f90
        fp16/lapack/cholqr/gecholqr_shift.f90

        fp16/lapack/householder/larf.f90
        fp16/lapack/householder/larfb.f90
        fp16/lapack/householder/larfbv.f90
        fp16/lapack/householder/larft.f90
        fp16/lapack/householder/larftv.f90
        fp16/lapack/householder/larf1f.f90
        fp16/lapack/householder/larfg.f90
        fp16/lapack/householder/larfgv1.f90
        fp16/lapack/householder/tprfb.f90
        fp16/lapack/qr/geqr2.f90
        fp16/lapack/qr/geqr2v.f90
        fp16/lapack/qr/geqrf.f90
        fp16/lapack/qr/geqrfv.f90
        fp16/lapack/qr/geqrt.f90
        fp16/lapack/qr/geqrt2.f90
        fp16/lapack/qr/geqrt3.f90
        fp16/lapack/qr/gemqrt.f90
        fp16/lapack/qr/geqrtv.f90
        fp16/lapack/qr/geqrt2v.f90
        fp16/lapack/qr/geqrt3v.f90
        fp16/lapack/qr/gemqrtv.f90
        fp16/lapack/qr/orm2r.f90
        fp16/lapack/qr/ormqr.f90
        fp16/lapack/qr/ormqrv.f90
        fp16/lapack/qr/orm2rv.f90
        fp16/lapack/qr/latsqr.f90
        fp16/lapack/qr/lamtsqr.f90
        fp16/lapack/qr/tpqrt.f90
        fp16/lapack/qr/tpqrt2.f90
        fp16/lapack/qr/tpmqrt.f90
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
        bf16/lpf_blas_bf16.F90)
endif()


SET(SRC ${SRC}
    highprecision/scale_diag_twoside.f90
    highprecision/scale_diag_right.f90
)

ADD_LIBRARY(lpf_blas
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level1>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level2>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_level3>>
    $<$<BOOL:${LPF_HAVE_FP16}>:$<TARGET_OBJECTS:fp16_aux>>

    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level1>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level2>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_level3>>
    $<$<BOOL:${LPF_HAVE_BF16}>:$<TARGET_OBJECTS:bf16_aux>>

    ${SRC}
)

TARGET_LINK_LIBRARIES(lpf_blas
    PRIVATE
    $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    PUBLIC
    LPF::LPF
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

ADD_LIBRARY(LPF::BLAS ALIAS lpf_blas)

INSTALL(TARGETS lpf_blas
    EXPORT  lpf_targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)



