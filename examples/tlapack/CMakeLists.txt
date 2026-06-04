include(FetchContent)
enable_language(CXX)
SET(BUILD_TESTING OFF)
FetchContent_Declare(
    tlapack
    GIT_REPOSITORY https://github.com/tlapack/tlapack.git
    GIT_TAG        v0.1.2
    EXCLUDE_FROM_ALL
)

FetchContent_MakeAvailable(tlapack)

add_executable(lu example_lu.cpp)
add_coverage(lu)
target_link_libraries(lu PRIVATE tlapack LPF::lpf LPF::lpf_blas)

add_executable(lu_fortran lu_call.f90 lu_interface.cpp)
add_coverage(lu_fortran)
target_link_libraries(lu_fortran PRIVATE tlapack LPF::lpf LPF::lpf_blas)

