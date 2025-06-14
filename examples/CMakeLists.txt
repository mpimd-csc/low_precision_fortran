SET (SRC
    demo1.f03
    demo_fraction.f03
)

FOREACH( S IN LISTS SRC)
    cmake_path(GET S STEM LAST_ONLY check_name)
    ADD_EXECUTABLE(${check_name} ${S})

    TARGET_LINK_LIBRARIES(${check_name}
        PRIVATE
            lpf
            $<$<BOOL:NEED_LINKING_AGAINST_LIBM>:${LIBM_NAME}>
    )

ENDFOREACH()


