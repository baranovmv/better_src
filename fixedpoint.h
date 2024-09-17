#ifndef SRC_FIXEDPOINT_H
#define SRC_FIXEDPOINT_H

#include <ostream>
#include <math.h>

template <typename T, unsigned int FractionalBits>
class FixedPoint
{
public:
    FixedPoint()
    : value(0)
    {}

    FixedPoint(T val)
    : value(val << FractionalBits)
    {}

    // Conversion from floating point to fixed point
    FixedPoint(float val)
    : value(static_cast<T>(val * (1 << FractionalBits)))
    {}

    template <unsigned int FractionalBitsN>
    FixedPoint(const FixedPoint &rhs)
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
        return static_cast<float>(value) / (1 << FractionalBits);
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

    FixedPoint operator-(const FixedPoint& other) const
    {
        return FixedPoint(value - other.value);
    }

    FixedPoint operator*(const FixedPoint& other) const
    {
        return FixedPoint((static_cast<T>(value) * other.value) >> FractionalBits);
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
    FixedPoint<T, FractionalBitsM> scaledown_resolution(float &residual)
    {
        static_assert(FractionalBitsM < FractionalBits);
        static constexpr size_t residual_bitlen = FractionalBits - FractionalBitsM;
        static constexpr size_t residual_mask = (1 << residual_bitlen) - 1;
        residual = static_cast<float>(value & residual_mask) / (1 << residual_bitlen);
        return FixedPoint<T, FractionalBitsM>(*this);
    }

    T get() const
    {
        return value;
    }

    /// Does linear interpolation between two values weighted by fractional part of value only.
    float fract_linear_interp(float x1, float x2)
    {
        const float fract = static_cast<float>(value & fract_bitmask_) / (1 << FractionalBits);
        return (x2 - x1) * fract + x1;
    }

private:
    T value;
    static constexpr T fract_bitmask_ = (1 << FractionalBits) - 1;
};

#endif //SRC_FIXEDPOINT_H
