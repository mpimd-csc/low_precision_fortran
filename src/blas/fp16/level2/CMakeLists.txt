SET(SRC
    hgbmv.c
    hgemv.c
    hger.c
    hsbmv.c
    hspmv.c
    hspr2.c
    hspr.c
    hsymv.c
    hsyr2.c
    hsyr.c
    htbmv.c
    htbsv.c
    htpmv.c
    htpsv.c
    htrmv.c
    htrsv.c

)

ADD_LIBRARY(fp16_level2 OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET fp16_level2 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(fp16_level2 PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)

TARGET_LINK_LIBRARIES(fp16_level2 PRIVATE LPF::LPF)

TARGET_INCLUDE_DIRECTORIES(fp16_level2
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src/fp16
)
