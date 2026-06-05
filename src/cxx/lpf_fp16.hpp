#ifndef LPF_FP16_HPP
#define LPF_FP16_HPP



#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <cmath>

struct lpf_fp16_t {
    _Float16 value;

    // Constructors
    lpf_fp16_t() : value(static_cast<_Float16>(0.0)) {}
    lpf_fp16_t(double val) : value(static_cast<_Float16>(val)) {}
    lpf_fp16_t(float val) : value(static_cast<_Float16>(val)) {}
    lpf_fp16_t(int val) : value(static_cast<_Float16>(val)) {}
    lpf_fp16_t(_Float16 val) : value(val) {}
    lpf_fp16_t(const lpf_fp16_t &val) : value(val.value) {}

    // Overloaded operators handling  lpf_fp16_t x lpf_fp16_t
    lpf_fp16_t operator+(const lpf_fp16_t& other) const {
        return lpf_fp16_t(value + other.value);
    }

    lpf_fp16_t operator-(const lpf_fp16_t& other) const {
        return lpf_fp16_t(value - other.value);
    }

    lpf_fp16_t operator*(const lpf_fp16_t& other) const {
        return lpf_fp16_t(value * other.value);
    }

    lpf_fp16_t operator/(const lpf_fp16_t& other) const {
        if (other.value == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp16_t(value / other.value);
    }

    lpf_fp16_t operator-() const {
        return lpf_fp16_t(-this->value);
    }

    // Overloaded operators handling  lpf_fp16_t x double
    lpf_fp16_t operator+(double other) const {
        return lpf_fp16_t(value + static_cast<_Float16>(other));
    }

    lpf_fp16_t operator-(double other) const {
        return lpf_fp16_t(value - static_cast<_Float16>(other));
    }

    lpf_fp16_t operator*(double other) const {
        return lpf_fp16_t(value * static_cast<_Float16>(other));
    }

    lpf_fp16_t operator/(double other) const {
        if (other == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp16_t(value / static_cast<_Float16>(other));
    }

    // Overloaded operators handling  lpf_fp16_t x float
    lpf_fp16_t operator+(float other) const {
        return lpf_fp16_t(value + static_cast<_Float16>(other));
    }

    lpf_fp16_t operator-(float other) const {
        return lpf_fp16_t(value - static_cast<_Float16>(other));
    }

    lpf_fp16_t operator*(float other) const {
        return lpf_fp16_t(value * static_cast<_Float16>(other));
    }

    lpf_fp16_t operator/(float other) const {
        if (other == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        return lpf_fp16_t(value / static_cast<_Float16>(other));
    }

    // Overloaded operators handling  lpf_fp16_t x= lpf_fp16_t
    lpf_fp16_t& operator+=(const lpf_fp16_t& other) {
        value += other.value;
        return *this;
    }

    lpf_fp16_t& operator-=(const lpf_fp16_t& other) {
        value -= other.value;
        return *this;
    }

    lpf_fp16_t& operator*=(const lpf_fp16_t& other) {
        value *= other.value;
        return *this;
    }

    lpf_fp16_t& operator/=(const lpf_fp16_t& other) {
        if (other.value == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        value /= other.value;
        return *this;
    }

    // Overloaded operators handling  lpf_fp16_t x= double
    lpf_fp16_t& operator+=(double other) {
        value += static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator-=(double other) {
        value -= static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator*=(double other) {
        value *= static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator/=(double other) {
        if (other == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        value /= static_cast<_Float16>(other);
        return *this;
    }

    // Overloaded operators handling  lpf_fp16_t x= float
    lpf_fp16_t& operator+=(float other) {
        value += static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator-=(float other) {
        value -= static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator*=(float other) {
        value *= static_cast<_Float16>(other);
        return *this;
    }

    lpf_fp16_t& operator/=(float other) {
        if (other == 0.0f) {
            throw std::runtime_error("Division by zero");
        }
        value /= static_cast<_Float16>(other);
        return *this;
    }

    // Overloaded Comparison operators lpf_fp16_t op lpf_fp16_t
    bool operator==(const lpf_fp16_t& other) const {
        return value == other.value;
    }

    bool operator!=(const lpf_fp16_t& other) const {
        return !(*this == other);
    }

    bool operator<(const lpf_fp16_t& other) const {
        return value < other.value;
    }

    bool operator<=(const lpf_fp16_t& other) const {
        return value <= other.value;
    }

    bool operator>(const lpf_fp16_t& other) const {
        return value > other.value;
    }

    bool operator>=(const lpf_fp16_t& other) const {
        return value >= other.value;
    }

    // Overloaded Comparison operators lpf_fp16_t op double
    bool operator==(double other) const {
        return value == static_cast<_Float16>(other);
    }

    bool operator!=(double other) const {
        return !(*this == other);
    }

    bool operator<(double other) const {
        return value < static_cast<_Float16>(other);
    }

    bool operator<=(double other) const {
        return value <= static_cast<_Float16>(other);
    }

    bool operator>(double other) const {
        return value > static_cast<_Float16>(other);
    }

    bool operator>=(double other) const {
        return value >= static_cast<_Float16>(other);
    }

    // Overloaded Comparison operators lpf_fp16_t op float
    bool operator==(float other) const {
        return value == static_cast<_Float16>(other);
    }

    bool operator!=(float other) const {
        return !(*this == other);
    }

    bool operator<(float other) const {
        return value < static_cast<_Float16>(other);
    }

    bool operator<=(float other) const {
        return value <= static_cast<_Float16>(other);
    }

    bool operator>(float other) const {
        return value > static_cast<_Float16>(other);
    }

    bool operator>=(float other) const {
        return value >= static_cast<_Float16>(other);
    }

    /*
       lpf_fp16_t& operator++() {  // Prefix
       ++value;
       return *this;
       }

       lpf_fp16_t operator++(int) {  // Postfix
       lpf_fp16_t temp(*this);
       ++value;
       return temp;
       }

       lpf_fp16_t& operator--() {  // Prefix
       --value;
       return *this;
       }

       lpf_fp16_t operator--(int) {  // Postfix
       lpf_fp16_t temp(*this);
       --value;
       return temp;
       }
       */

    // Overloaded IO
    friend std::ostream& operator<<(std::ostream& os, const lpf_fp16_t& f) {
        double tmp = static_cast<double>(f.value);
        os << tmp;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, lpf_fp16_t& f) {
        double tmp;
        is >> tmp;
        f = tmp;
        return is;
    }

    // Typecasts
    operator double() const {
        return static_cast<double>(value);
    }

    // Typumwandlung in float
    operator float() const {
        return static_cast<float>(value);
    }

    // Typumwandlung in int
    operator int() const {
        return static_cast<int>(value);
    }
};

// Arithmetics for left sided operators.
inline lpf_fp16_t operator+(double left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) + right.value);
}

inline lpf_fp16_t operator-(double left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) - right.value);
}

inline lpf_fp16_t operator*(double left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) * right.value);
}

inline lpf_fp16_t operator/(double left, const lpf_fp16_t& right) {
    if (right.value == 0.0) {
        throw std::runtime_error("Division by zero");
    }
    return lpf_fp16_t(static_cast<_Float16>(left) / right.value);
}

inline lpf_fp16_t operator+(float left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) + right.value);
}

inline lpf_fp16_t operator-(float left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) - right.value);
}

inline lpf_fp16_t operator*(float left, const lpf_fp16_t& right) {
    return lpf_fp16_t(static_cast<_Float16>(left) * right.value);
}

inline lpf_fp16_t operator/(float left, const lpf_fp16_t& right) {
    if (right.value == 0.0) {
        throw std::runtime_error("Division by zero");
    }
    return lpf_fp16_t(static_cast<_Float16>(left) / right.value);
}

// Comparison Operators
inline bool operator==(double left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) == right.value;
}

inline bool operator!=(double left, const lpf_fp16_t& right) {
    return !(left == right);
}

inline bool operator<(double left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) < right.value;
}

inline bool operator<=(double left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) <= right.value;
}

inline bool operator>(double left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) > right.value;
}

inline bool operator>=(double left, const lpf_fp16_t& right) {
    return left >= right.value;
}

inline bool operator==(float left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) == right.value;
}

inline bool operator!=(float left, const lpf_fp16_t& right) {
    return !(left == right);
}

inline bool operator<(float left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) < right.value;
}

inline bool operator<=(float left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) <= right.value;
}

inline bool operator>(float left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) > right.value;
}

inline bool operator>=(float left, const lpf_fp16_t& right) {
    return static_cast<_Float16>(left) >= right.value;
}

namespace std {
    inline const lpf_fp16_t abs   (const lpf_fp16_t& x) {
        if (static_cast<double>(x.value) >= 0.0) return x;
        else return -x;
    }

    inline bool isinf(const lpf_fp16_t& x){
        return std::isinf(static_cast<double>(x.value));
    }

    inline bool isnan(const lpf_fp16_t& x){
        return std::isnan(static_cast<double>(x.value));
    }

    inline const lpf_fp16_t sqrt (const lpf_fp16_t& x) {
        return std::sqrt(static_cast<double>(x.value));
    }

    inline const lpf_fp16_t pow(const lpf_fp16_t& a, const lpf_fp16_t& b ){
        return lpf_fp16_t(std::pow(static_cast<double>(a.value), static_cast<double>(b.value)));
    }
    inline const lpf_fp16_t pow(const int a, const lpf_fp16_t& b ){
        return lpf_fp16_t(std::pow(static_cast<double>(a), static_cast<double>(b.value)));
    }

    inline const lpf_fp16_t log2(const lpf_fp16_t& x) {
        return std::log2(static_cast<double>(x.value));
    }

    inline const lpf_fp16_t ceil(const lpf_fp16_t& x) {
        return std::ceil(static_cast<double>(x.value));
    }
    inline const lpf_fp16_t floor(const lpf_fp16_t& x) {
        return std::floor(static_cast<double>(x.value));
    }
}

namespace std
{

    template<>
        class numeric_limits<lpf_fp16_t>
        {
            public:
                static const bool is_specialized    = true;
                static const bool is_signed         = true;
                static const bool is_integer        = false;
                static const bool is_exact          = false;
                static const int  radix             = 2;
                static const int  min_exponent      = __FLT16_MIN_EXP__;
                static const int  max_exponent      = __FLT16_MAX_EXP__;

        };
}

#endif

