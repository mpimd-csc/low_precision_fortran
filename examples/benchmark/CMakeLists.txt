ENABLE_LANGUAGE(CXX)
set(CMAKE_CXX_STANDARD 17)
include(FetchContent)
FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG v4.3.11
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nanobench)

ADD_EXECUTABLE(run_benchmark linspace.f90 benchmark.cpp)
add_coverage(run_benchmark)
target_link_libraries(run_benchmark PRIVATE LPF::lpf nanobench)

