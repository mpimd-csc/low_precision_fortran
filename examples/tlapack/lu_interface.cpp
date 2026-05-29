
#include "lpf_fp8_e4m3.hpp"
#include "tlapack/base/types.hpp"
#include <tlapack/plugins/legacyArray.hpp>
#include <lpf_fp16_tlapack.hpp>
#include <lpf_bf16_tlapack.hpp>
#include <lpf_fp8_e4m3_tlapack.hpp>
#include <lpf_fp8_e5m2_tlapack.hpp>

// <T>LAPACK
#include <tlapack/blas/trsv.hpp>
#include <tlapack/lapack/getrf.hpp>
#include <tlapack/lapack/lacpy.hpp>
#include <tlapack/lapack/lange.hpp>

// C++ headers
#include <cstdint>
#include <iostream>

    template<typename T>
void solve (int64_t n, T * A, int64_t lda, T *x, int64_t incx)
{
    const tlapack::Layout L = tlapack::Layout::ColMajor;
    using real_t = tlapack::real_type<T>;
    using idx_t = size_t;

    // Create the n-by-n matrix A
    tlapack::LegacyMatrix<T, idx_t, L> tA(n, n, A, lda);
    tlapack::LegacyVector<T, idx_t> tx(n, x, incx);

    // Allocate space for the LU decomposition
    std::vector<size_t> piv(n);
    std::vector<T> LU_(n * n);
    tlapack::LegacyMatrix<T, idx_t, L> LU(n, n, LU_.data(), n);

    // Matrix A is kept unchanged
    tlapack::lacpy(tlapack::GENERAL, tA, LU);

    // Computing the LU decomposition of A
    int info = tlapack::getrf(LU, piv);
    if (info != 0) {
        std::cerr << "Matrix could not be factorized!" << std::endl;
        return;
    }

    for (idx_t i = 0; i < n ; i++) {
        if (piv[i] != i) {
            auto val1 = tx[i];
            tx[i] = tx[piv[i]];
            tx[piv[i]] = val1;
        }
    }
    tlapack::trsv(tlapack::Uplo::Lower, tlapack::Op::NoTrans, tlapack::Diag::Unit, LU, tx);
    tlapack::trsv(tlapack::Uplo::Upper, tlapack::Op::NoTrans, tlapack::Diag::NonUnit, LU, tx);
}

    extern "C"
void solve_fp16(int64_t n, lpf_fp16_t * A, int64_t lda, lpf_fp16_t *x, int64_t incx)
{
    solve<lpf_fp16_t>(n, A, lda, x, incx);
}

    extern "C"
void solve_bf16(int64_t n, lpf_bf16_t * A, int64_t lda, lpf_bf16_t *x, int64_t incx)
{
    solve<lpf_bf16_t>(n, A, lda, x, incx);
}

    extern "C"
void solve_fp8_e4m3(int64_t n, lpf_fp8_e4m3_t * A, int64_t lda, lpf_fp8_e4m3_t *x, int64_t incx)
{
    solve<lpf_fp8_e4m3_t>(n, A, lda, x, incx);
}
    extern "C"
void solve_fp8_e5m2(int64_t n, lpf_fp8_e5m2_t * A, int64_t lda, lpf_fp8_e5m2_t *x, int64_t incx)
{
    solve<lpf_fp8_e5m2_t>(n, A, lda, x, incx);
}

