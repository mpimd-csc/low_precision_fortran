# Examples

This directory contains examples demonstrating the use of the low-precision types provided by the LPF library.

## Basic Type Demos

These standalone Fortran files provide a quick way to verify the basic functionality and properties of the supported low-precision types (`BF16`, `FP16`, and `FP8`).

- **Property Demos**: `bf16_demo1.f90` and `fp16_demo1.f90` demonstrate basic properties such as `EPSILON`, `HUGE`, `TINY`, and exponent ranges.
- **Component Extraction**: `bf16_demo_fraction.f90` and `fp16_demo_fraction.f90` show how to extract the exponent and fraction components from a value.
- **Formatting and Stats**: `bf16_demo_stats.f90` and `fp16_demo_stats.f90` demonstrate the use of the `DT` output format and basic functions like `ABS`.
- **Simple BLAS**: `fp8_blas.f90` is a simple example demonstrating the use of `asum` for the `FP8_E5M2` type.

## Algorithmic Examples

- **Conjugate Gradient (`blas_cg.f90`)**: This example implements the Conjugate Gradient (CG) iterative solver. It demonstrates how to use LPF's BLAS-like operations (`gemv`, `dot`, `axpy`, `nrm2`) to solve linear systems across multiple precisions (`BF16`, `FP16`, `FP8_E4M3`, and `FP8_E5M2`).

## Performance Benchmarking

The `benchmark/` directory contains micro-benchmarks to evaluate the performance of low-precision operations. It uses a C++ driver with the `nanobench` library to measure the execution time and performance counters of Fortran routines (e.g., `linspace.f90`).

## <T>LAPACK Integration

The `tlapack/` directory demonstrates how LPF types can be integrated with the $\langle T \rangle$LAPACK library for high-level linear algebra operations. For example, `example_lu.cpp` and `lu_call.f90` show how to compute the inverse of a matrix via LU decomposition using `FP16`, `BF16`, and `FP8` precisions.
