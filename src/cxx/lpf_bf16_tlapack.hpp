#ifndef LPF_BF16_TLAPACK_HPP
#define LPF_BF16_TLAPACK_HPP

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "lpf_bf16.hpp"

#include "tlapack/base/types.hpp"

namespace tlapack {

    namespace traits {
        // mpfr::mpreal is a real type that satisfies tlapack::concepts::Real
        template <>
            struct real_type_traits<lpf_bf16_t, int> {
                using type = lpf_bf16_t;
                constexpr static bool is_real = true;
            };
        // The complex type of mpfr::mpreal is std::complex<mpfr::mpreal>
        template <>
            struct complex_type_traits<lpf_bf16_t, int> {
                using type = std::complex<lpf_bf16_t>;
                constexpr static bool is_complex = false;
            };
    }  // namespace traits

    // Argument-dependent lookup (ADL) will include the remaining functions,
    // e.g., mpfr::sin, mpfr::cos.
    // Including them here may cause ambiguous call of overloaded function.
    // See: https://en.cppreference.com/w/cpp/language/adl

    // Forward declaration
    template <typename real_t>
        int digits() noexcept;

    // Specialization for the mpfr::mpreal datatype
    template <>
        inline int digits<lpf_bf16_t>() noexcept
        {
            // return std::numeric_limits<mpfr::mpreal>::digits();
            return 10;
        }


}  // namespace tlapack


#endif
