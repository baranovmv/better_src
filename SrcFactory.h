#pragma once
#include <memory>

#include "converter.h"

// Define the enum with 5 values
enum class SAMPLE_TYPE {
    S8,
    S16,
    S20,
    S24,
    F32,
};

namespace factory_detail
{
    // Helper template for generating the correct class based on parameters
    template<int n_channels, typename... Args>
    ISrc* createImpl(SAMPLE_TYPE sample_type, Args&&... args) {
        switch (sample_type) {
        case SAMPLE_TYPE::F32:
            return new Src<n_channels, float, float>(std::forward<Args>(args)...);
        case SAMPLE_TYPE::S16:
            return new Src<n_channels, FixedPoint<int32_t, int64_t, 16>, FixedPoint<int64_t, int64_t, 32>>(std::forward<Args>(args)...);
        default:
            return nullptr;
        }
    }

    // Factory implementation that uses template recursion
    template<int n_channels, typename = void>
    struct FactoryImpl {
        template<typename... Args>
        static ISrc* create(int n, SAMPLE_TYPE sample_width, Args&&... args)
        {
            if (n == n_channels) {
                return createImpl<n_channels>(sample_width, std::forward<Args>(args)...);
            }
            return FactoryImpl<n_channels-1>::create(n, sample_width, std::forward<Args>(args)...);
        }
    };

    // Base case specialization
    template<typename Dummy>
    struct FactoryImpl<0, Dummy> {
        template<typename... Args>
        static ISrc* create(int n_channels, SAMPLE_TYPE sample_width, Args&&... args) {
            throw std::invalid_argument("Invalid parameters");
        }
    };
}

class SrcFactory {
public:
    template<typename... Args>
    static ISrc* create(int n_channels, SAMPLE_TYPE sample_width, Args&&... args) {
        if (n_channels <= 0 || n_channels > 2) {
            throw std::invalid_argument("n_channels must be between 1 and 2");
        }
        return factory_detail::FactoryImpl<2>::create(n_channels, sample_width, std::forward<Args>(args)...);
    }
};
