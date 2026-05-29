SET(SRC
    bgbmv.c
    bgemv.c
    bger.c
    bsbmv.c
    bspmv.c
    bspr2.c
    bspr.c
    bsymv.c
    bsyr2.c
    bsyr.c
    btbmv.c
    btbsv.c
    btpmv.c
    btpsv.c
    btrmv.c
    btrsv.c

)

ADD_LIBRARY(bf16_level2 OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET bf16_level2 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(bf16_level2 PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)

TARGET_LINK_LIBRARIES(bf16_level2 PRIVATE LPF::lpf)

TARGET_INCLUDE_DIRECTORIES(bf16_level2
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src/bf16
)
