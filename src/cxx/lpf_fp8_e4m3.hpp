#ifndef LPF_FP8_E4M3_HPP
#define LPF_FP8_E4M3_HPP



#include <iostream>
// #include <iomanip>
#include <stdexcept>
#include <limits>
#include <cmath>
#include "lpf_fp8_e4m3_emu.h"

struct lpf_fp8_e4m3_t {
    fp8_e4m3_t value;

    // Constructors
    lpf_fp8_e4m3_t() {
        this->value = fp8_e4m3_from_float(0.0f);
    }
    lpf_fp8_e4m3_t(double val) {
        this->value = fp8_e4m3_from_float(static_cast<float>(val));
    }
    lpf_fp8_e4m3_t(float val) {
        this->value = fp8_e4m3_from_float(val);
    }
    lpf_fp8_e4m3_t(int val) {
        this->value = fp8_e4m3_from_float(static_cast<float>(val));
    }
    lpf_fp8_e4m3_t(fp8_e4m3_t val) {
        this->value = val;
    }
    lpf_fp8_e4m3_t(const lpf_fp8_e4m3_t &val) {
        this->value = val.value;
    }

    // Overloaded operators handling  lpf_fp8_e4m3_t x lpf_fp8_e4m3_t
    lpf_fp8_e4m3_t operator+(const lpf_fp8_e4m3_t& other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_add(value, other.value));
    }

    lpf_fp8_e4m3_t operator-(const lpf_fp8_e4m3_t& other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_sub(value,other.value));
    }

    lpf_fp8_e4m3_t operator*(const lpf_fp8_e4m3_t& other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_mul(value, other.value));
    }

    lpf_fp8_e4m3_t operator/(const lpf_fp8_e4m3_t& other) const {
        if (fp8_e4m3_to_float(other.value) == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp8_e4m3_t(fp8_e4m3_div(value, other.value));
    }

    lpf_fp8_e4m3_t operator-() const {
        return lpf_fp8_e4m3_t(fp8_e4m3_unitary_minus(this->value));
    }

    // Overloaded operators handling  lpf_fp8_e4m3_t x double
    lpf_fp8_e4m3_t operator+(const double other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_add(value, fp8_e4m3_from_float(static_cast<float>(other))));
    }

    lpf_fp8_e4m3_t operator-(const double other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_sub(value, fp8_e4m3_from_float(static_cast<float>(other))));
    }

    lpf_fp8_e4m3_t operator*(const double other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_mul(value, fp8_e4m3_from_float(static_cast<float>(other))));
    }

    lpf_fp8_e4m3_t operator/(const double other) const {
        if (other == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp8_e4m3_t(fp8_e4m3_div(value, fp8_e4m3_from_float(static_cast<float>(other))));
    }



    // Overloaded operators handling  lpf_fp8_e4m3_t x float
    lpf_fp8_e4m3_t operator+(const float other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_add(value, fp8_e4m3_from_float(other)));
    }

    lpf_fp8_e4m3_t operator-(const float other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_sub(value, fp8_e4m3_from_float(other)));
    }

    lpf_fp8_e4m3_t operator*(const float other) const {
        return lpf_fp8_e4m3_t(fp8_e4m3_mul(value, fp8_e4m3_from_float(other)));
    }

    lpf_fp8_e4m3_t operator/(const float other) const {
        if (other == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp8_e4m3_t(fp8_e4m3_div(value, fp8_e4m3_from_float(other)));
    }



    // Overloaded operators handling  lpf_fp8_e4m3_t x= lpf_fp8_e4m3_t
    lpf_fp8_e4m3_t& operator+=(const lpf_fp8_e4m3_t& other) {
        value = fp8_e4m3_add(value, other.value);
        return *this;
    }

    lpf_fp8_e4m3_t& operator-=(const lpf_fp8_e4m3_t& other) {
        value = fp8_e4m3_sub(value, other.value);
        return *this;
    }

    lpf_fp8_e4m3_t& operator*=(const lpf_fp8_e4m3_t& other) {
        value = fp8_e4m3_mul(value, other.value);
        return *this;
    }

    lpf_fp8_e4m3_t& operator/=(const lpf_fp8_e4m3_t& other) {
        if (fp8_e4m3_to_float(other.value) == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        value = fp8_e4m3_div(value, other.value);
        return *this;
    }

    // Overloaded operators handling  lpf_fp8_e4m3_t x= double
    lpf_fp8_e4m3_t& operator+=(double other) {
        value = fp8_e4m3_add(value, fp8_e4m3_from_float(static_cast<float>(other)));
        return *this;
    }

    lpf_fp8_e4m3_t& operator-=(double other) {
        value = fp8_e4m3_sub(value, fp8_e4m3_from_float(static_cast<float>(other)));
        return *this;
    }

    lpf_fp8_e4m3_t& operator*=(double other) {
        value = fp8_e4m3_mul(value, fp8_e4m3_from_float(static_cast<float>(other)));
        return *this;
    }

    lpf_fp8_e4m3_t& operator/=(double other) {
        if (other == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        value = fp8_e4m3_div(value, fp8_e4m3_from_float(static_cast<float>(other)));
        return *this;
    }

    // Overloaded operators handling  lpf_fp8_e4m3_t x= float
    lpf_fp8_e4m3_t& operator+=(float other) {
        value = fp8_e4m3_add(value, fp8_e4m3_from_float(other));
        return *this;
    }

    lpf_fp8_e4m3_t& operator-=(float other) {
        value = fp8_e4m3_sub(value, fp8_e4m3_from_float(other));
        return *this;
    }

    lpf_fp8_e4m3_t& operator*=(float other) {
        value = fp8_e4m3_mul(value, fp8_e4m3_from_float(other));
        return *this;
    }

    lpf_fp8_e4m3_t& operator/=(float other) {
        if (other == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        value = fp8_e4m3_div(value, fp8_e4m3_from_float(other));
        return *this;
    }

    // Overloaded Comparison operators lpf_fp8_e4m3_t op lpf_fp8_e4m3_t
    bool operator==(const lpf_fp8_e4m3_t& other) const {
        return value == other.value;
    }

    bool operator!=(const lpf_fp8_e4m3_t& other) const {
        return !(*this == other);
    }

    bool operator<(const lpf_fp8_e4m3_t& other) const {
        return fp8_e4m3_to_float(value) < fp8_e4m3_to_float(other.value);
    }

    bool operator<=(const lpf_fp8_e4m3_t& other) const {
        return fp8_e4m3_to_float(value) <= fp8_e4m3_to_float(other.value);
    }

    bool operator>(const lpf_fp8_e4m3_t& other) const {
        return fp8_e4m3_to_float(value) > fp8_e4m3_to_float(other.value);
    }

    bool operator>=(const lpf_fp8_e4m3_t& other) const {
        return fp8_e4m3_to_float(value) >= fp8_e4m3_to_float(other.value);
    }

    // Overloaded Comparison operators lpf_fp8_e4m3_t op double
    bool operator==(double other) const {
        return fp8_e4m3_to_float(value) == static_cast<float>(other);
    }

    bool operator!=(double other) const {
        return !(*this == other);
    }

    bool operator<(double other) const {
        return fp8_e4m3_to_float(value) < static_cast<float>(other);
    }

    bool operator<=(double other) const {
        return fp8_e4m3_to_float(value) <= static_cast<float>(other);
    }

    bool operator>(double other) const {
        return fp8_e4m3_to_float(value) > static_cast<float>(other);
    }

    bool operator>=(double other) const {
        return fp8_e4m3_to_float(value) >= static_cast<float>(other);
    }

    // Overloaded Comparison operators lpf_fp8_e4m3_t op float
    bool operator==(float other) const {
        return fp8_e4m3_to_float(value) == other;
    }

    bool operator!=(float other) const {
        return !(*this == other);
    }

    bool operator<(float other) const {
        return fp8_e4m3_to_float(value) < other;
    }

    bool operator<=(float other) const {
        return fp8_e4m3_to_float(value) <= other;
    }

    bool operator>(float other) const {
        return fp8_e4m3_to_float(value) > other;
    }

    bool operator>=(float other) const {
        return fp8_e4m3_to_float(value) >= other;
    }

    /*
       lpf_fp8_e4m3_t& operator++() {  // Prefix
       ++value;
       return *this;
       }

       lpf_fp8_e4m3_t operator++(int) {  // Postfix
       lpf_fp8_e4m3_t temp(*this);
       ++value;
       return temp;
       }

       lpf_fp8_e4m3_t& operator--() {  // Prefix
       --value;
       return *this;
       }

       lpf_fp8_e4m3_t operator--(int) {  // Postfix
       lpf_fp8_e4m3_t temp(*this);
       --value;
       return temp;
       }
       */

    // Overloaded IO
    friend std::ostream& operator<<(std::ostream& os, const lpf_fp8_e4m3_t& f) {
        float tmp = fp8_e4m3_to_float(f.value);
        os << tmp;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, lpf_fp8_e4m3_t& f) {
        float tmp;
        is >> tmp;
        f.value  = fp8_e4m3_from_float(tmp);
        return is;
    }

    // Typecasts
    operator double() const {
        return static_cast<double>(fp8_e4m3_to_float(value));
    }

    // Typumwandlung in float
    operator float() const {
        return static_cast<float>(fp8_e4m3_to_float(value));
    }

    // Typumwandlung in int
    operator int() const {
        return static_cast<int>(fp8_e4m3_to_float(value));
    }
};

// Arithmetics for left sided operators.
inline lpf_fp8_e4m3_t operator+(double left, const lpf_fp8_e4m3_t& right) {
    return lpf_fp8_e4m3_t(fp8_e4m3_add(fp8_e4m3_from_float(static_cast<float>(left)),right.value ));
}

inline lpf_fp8_e4m3_t operator-(double left, const lpf_fp8_e4m3_t& right) {
    return lpf_fp8_e4m3_t(fp8_e4m3_sub(fp8_e4m3_from_float(static_cast<float>(left)),right.value ));
}

inline lpf_fp8_e4m3_t operator*(double left, const lpf_fp8_e4m3_t& right) {
    return lpf_fp8_e4m3_t(fp8_e4m3_mul(fp8_e4m3_from_float(static_cast<float>(left)),right.value ));
}

inline lpf_fp8_e4m3_t operator/(double left, const lpf_fp8_e4m3_t& right) {
    if (fp8_e4m3_to_float(right.value) == 0.0f) {
        throw std::runtime_error("Division by zero");
    }
    return lpf_fp8_e4m3_t(fp8_e4m3_div(fp8_e4m3_from_float(static_cast<float>(left)),right.value ));
}

inline lpf_fp8_e4m3_t operator+(float left, const lpf_fp8_e4m3_t& right) {
    return lpf_fp8_e4m3_t(fp8_e4m3_add(fp8_e4m3_from_float(left),right.value ));
}

inline lpf_fp8_e4m3_t operator-(float left, const lpf_fp8_e4m3_t& right) {
    return lpf_fp8_e4m3_t(fp8_e4m3_sub(fp8_e4m3_from_float(left),right.value ));
}

inline lpf_fp8_e4m3_t operator*(float left, const lpf_fp8_e4m3_t& right) {

    return lpf_fp8_e4m3_t(fp8_e4m3_mul(fp8_e4m3_from_float(left),right.value ));
}

inline lpf_fp8_e4m3_t operator/(float left, const lpf_fp8_e4m3_t& right) {
    if (fp8_e4m3_to_float(right.value) == 0.0f) {
        throw std::runtime_error("Division by zero");
    }
    return lpf_fp8_e4m3_t(fp8_e4m3_div(fp8_e4m3_from_float(left),right.value ));
}

// Comparison Operators
inline bool operator==(double left, const lpf_fp8_e4m3_t& right) {
    return static_cast<float>(left) == fp8_e4m3_to_float(right.value);
}

inline bool operator!=(double left, const lpf_fp8_e4m3_t& right) {
    return !(left == right);
}

inline bool operator<(double left, const lpf_fp8_e4m3_t& right) {

    return static_cast<float>(left) < fp8_e4m3_to_float(right.value);
}

inline bool operator<=(double left, const lpf_fp8_e4m3_t& right) {

    return static_cast<float>(left) <= fp8_e4m3_to_float(right.value);
}

inline bool operator>(double left, const lpf_fp8_e4m3_t& right) {

    return static_cast<float>(left) > fp8_e4m3_to_float(right.value);
}

inline bool operator>=(double left, const lpf_fp8_e4m3_t& right) {

    return static_cast<float>(left) >= fp8_e4m3_to_float(right.value);
}

inline bool operator==(float left, const lpf_fp8_e4m3_t& right) {
    return left == fp8_e4m3_to_float(right.value);
}

inline bool operator!=(float left, const lpf_fp8_e4m3_t& right) {
    return !(left == right);
}

inline bool operator<(float left, const lpf_fp8_e4m3_t& right) {
    return left < fp8_e4m3_to_float(right.value);
}

inline bool operator<=(float left, const lpf_fp8_e4m3_t& right) {
    return left <= fp8_e4m3_to_float(right.value);
}

inline bool operator>(float left, const lpf_fp8_e4m3_t& right) {
    return left > fp8_e4m3_to_float(right.value);
}

inline bool operator>=(float left, const lpf_fp8_e4m3_t& right) {
    return left >= fp8_e4m3_to_float(right.value);
}

namespace std {
    inline const lpf_fp8_e4m3_t abs   (const lpf_fp8_e4m3_t& x) {
        return lpf_fp8_e4m3_t(fp8_e4m3_abs(x.value));
    }

    inline bool isinf(const lpf_fp8_e4m3_t& x){
        return fp8_e4m3_isinf(x.value);
    }

    inline bool isnan(const lpf_fp8_e4m3_t& x){
        return fp8_e4m3_isnan(x.value);
    }

    inline const lpf_fp8_e4m3_t sqrt (const lpf_fp8_e4m3_t& x) {
        return lpf_fp8_e4m3_t(fp8_e4m3_sqrt(x.value));
    }

    inline const lpf_fp8_e4m3_t pow(const lpf_fp8_e4m3_t& a, const lpf_fp8_e4m3_t& b ){
        return lpf_fp8_e4m3_t(fp8_e4m3_pow(a.value, b.value));
    }
    inline const lpf_fp8_e4m3_t pow(const int a, const lpf_fp8_e4m3_t& b ){
        float _a = static_cast<float>(a);
        return lpf_fp8_e4m3_t(fp8_e4m3_pow(fp8_e4m3_from_float(_a), b.value));
    }

    inline const lpf_fp8_e4m3_t log2(const lpf_fp8_e4m3_t& x) {
        return lpf_fp8_e4m3_t(fp8_e4m3_log2(x.value));
    }

    inline const lpf_fp8_e4m3_t ceil(const lpf_fp8_e4m3_t& x) {
        float _x = fp8_e4m3_from_float(x);
        return lpf_fp8_e4m3_t(std::ceil(_x));
    }
    inline const lpf_fp8_e4m3_t floor(const lpf_fp8_e4m3_t& x) {
        float _x = fp8_e4m3_from_float(x);
        return lpf_fp8_e4m3_t(std::ceil(_x));
    }
}

namespace std
{

    template<>
        class numeric_limits<lpf_fp8_e4m3_t>
        {
            public:
                static const bool is_specialized    = true;
                static const bool is_signed         = true;
                static const bool is_integer        = false;
                static const bool is_exact          = false;
                static const int  radix             = 2;
                static const int  min_exponent      = __FP8_E4M3_MIN_EXP__;
                static const int  max_exponent      = __FP8_E4M3_MAX_EXP__;

        };
}

#endif

