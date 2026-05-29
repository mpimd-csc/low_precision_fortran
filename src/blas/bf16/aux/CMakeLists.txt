SET(SRC
    b2scopy.c
    b2slacpy.c
    blacpy.c
    blaset.c
    s2bcopy.c
    s2blacpy.c
)
ADD_LIBRARY(bf16_aux OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET bf16_aux PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(bf16_aux PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)

TARGET_LINK_LIBRARIES(bf16_aux PRIVATE LPF::lpf)

TARGET_INCLUDE_DIRECTORIES(bf16_aux
    PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/bf16
)
