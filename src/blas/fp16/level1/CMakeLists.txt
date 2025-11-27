SET(SRC
    hasum.c
    haxpby.c
    haxpy.c
    hcopy.c
    hdot.c
    hnrm2.c
    hnrm2_fp32.c
    hrot.c
    hrotg.c
    hrotm.c
    hrotmg.c
    hscal.c
    hshdot.c
    hswap.c
    ihamax.c
)

ADD_LIBRARY(fp16_level1 OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET fp16_level1 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(fp16_level1 PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)
TARGET_LINK_LIBRARIES(fp16_level1 PUBLIC LPF::LPF)

TARGET_INCLUDE_DIRECTORIES(fp16_level1
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src/fp16
)
