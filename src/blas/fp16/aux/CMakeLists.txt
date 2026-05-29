SET(SRC
    h2scopy.c
    h2slacpy.c
    hlacpy.c
    hlaset.c
    s2hcopy.c
    s2hlacpy.c
)
ADD_LIBRARY(fp16_aux OBJECT ${SRC})

if (LPF_USE_LTO)
    set_property(TARGET fp16_aux PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
TARGET_COMPILE_DEFINITIONS(fp16_aux PRIVATE
    $<$<BOOL:${LPF_INTEGER8}>:-DLPF_INTEGER8>
)

TARGET_LINK_LIBRARIES(fp16_aux PRIVATE LPF::lpf)

TARGET_INCLUDE_DIRECTORIES(fp16_aux
    PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/fp16
)
