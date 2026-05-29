SET(SRC
    basum.c
    baxpby.c
    baxpy.c
    bcopy.c
    bdot.c
    bnrm2.c
    bnrm2_fp32.c
    brot.c
    brotg.c
    brotm.c
    brotmg.c
    bscal.c
    bsbdot.c
    bswap.c
    ibamax.c
)

ADD_LIBRARY(bf16_level1 OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET bf16_level1 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(bf16_level1 PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)
TARGET_LINK_LIBRARIES(bf16_level1 PUBLIC LPF::lpf)

TARGET_INCLUDE_DIRECTORIES(bf16_level1
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src/bf16
)
