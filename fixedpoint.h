#ifndef SRC_FIXEDPOINT_H
#define SRC_FIXEDPOINT_H

#include <ostream>
#include <math.h>

template <typename T, typename LONG_T, unsigned int FractionalBits>
class FixedPoint
{
public:
    FixedPoint()
    : value(0)
    {}

    explicit FixedPoint(T val)
    : value(val << FractionalBits)
    {}

    explicit FixedPoint(LONG_T val)
    : value(val >> FractionalBits)
    {}

    // Conversion from floating point to fixed point
    explicit FixedPoint(float val)
    : value(static_cast<T>(val * (T(1) << FractionalBits)))
    {}

    template <unsigned int FractionalBitsN>
    explicit FixedPoint(const FixedPoint &rhs)
    {
        if (FractionalBitsN > FractionalBits) {
            value = rhs.value >> (FractionalBitsN - FractionalBits);
        } else {
            value = rhs.value << (FractionalBits - FractionalBitsN);
        }
    }

    // Conversion from fixed point to floating point
    operator float() const
    {
        return static_cast<float>(value) / (T(1) << FractionalBits);
    }

    inline T floor() const
    {
        return value >> FractionalBits;
    }

    inline T ceil() const
    {
        if ((value & fract_bitmask_) != 0) {
            return floor() + 1;
        } else {
            return floor();
        }
    }

    FixedPoint operator+(const FixedPoint& other) const
    {
        return FixedPoint(value + other.value);
    }

    FixedPoint operator+=(const FixedPoint& other)
    {
        value = value + other.value;
        return *this;
    }

    FixedPoint operator-(const FixedPoint& other) const
    {
        return FixedPoint(value - other.value);
    }

    FixedPoint operator-=(const FixedPoint& other)
    {
        value = value - other.value;
        return *this;
    }

    FixedPoint operator*(const FixedPoint& other) const
    {
        const LONG_T accum = LONG_T(value) * LONG_T(other.value);
        const auto result = FixedPoint(accum);
        return result;
    }

    FixedPoint operator/(const FixedPoint& other) const
    {
        return FixedPoint((static_cast<T>(value) << FractionalBits) / other.value);
    }

    friend std::ostream& operator<<(std::ostream& os, const FixedPoint& fp)
    {
        os << static_cast<float>(fp);
        return os;
    }

    template<unsigned int FractionalBitsM>
    FixedPoint<T, LONG_T, FractionalBitsM> scaledown_resolution(float &residual)
    {
        static_assert(FractionalBitsM < FractionalBits);
        static constexpr size_t residual_bitlen = FractionalBits - FractionalBitsM;
        static constexpr size_t residual_mask = (T(1) << residual_bitlen) - 1;
        residual = static_cast<float>(value & residual_mask) / (T(1) << residual_bitlen);
        return FixedPoint<T, LONG_T, FractionalBitsM>(*this);
    }

    T get() const
    {
        return value;
    }

    static FixedPoint<T, LONG_T, FractionalBits> convert(const T& x)
    {
        return FixedPoint<T, LONG_T, FractionalBits>(x);
    }

    /// Does linear interpolation between two values weighted by fractional part of value only.
    float fract_linear_interp(float x1, float x2)
    {
        const float fract = static_cast<float>(value & fract_bitmask_) / (T(1) << FractionalBits);
        return (x2 - x1) * fract + x1;
    }

private:
    T value;
    static constexpr T fract_bitmask_ = (T(1) << FractionalBits) - 1;
};

#endif //SRC_FIXEDPOINT_H
