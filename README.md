# HPBlas

A **H**alf-**P**recision **BLAS** implementation

## Installation

The library currently only works with gcc compilers supporting the `_Float16`
datatype and the `-m16fc` switch. This starts with GCC 11. For properhardware
support, a CPU with `AVX512FP16` is required.

### Building
```shell
cmake -S . -B build-dir -DCMAKE_INSTALL_PREFIX=INSTALL-LOCATION
cmake --build build-dir
cmake --build build-dir --target install
```

### Testing
```shell
(cd build-dir; ctest --output-on-failure)
```


