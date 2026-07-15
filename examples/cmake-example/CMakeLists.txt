cmake_minimum_required(VERSION 3.18)
project(cmake-example LANGUAGES C Fortran)

find_package(LPF REQUIRED)

add_executable(blas_example fp8_blas.f90)
target_link_libraries(blas_example PRIVATE LPF::lpf LPF::lpf_blas )
