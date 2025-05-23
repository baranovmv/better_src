#pragma once
#include <memory>

#include "converter.h"

namespace factory_detail
{
    // Helper template for generating the correct class based on parameters
    template<int n_channels, typename... Args>
    ISrc* createImpl(Args&&... args) {
        return new Src<n_channels>(std::forward<Args>(args)...);
    }

    // Factory implementation that uses template recursion
    template<int n_channels, typename = void>
    struct FactoryImpl {
        template<typename... Args>
        static ISrc* create(int n, Args&&... args)
        {
            if (n == n_channels) {
                return createImpl<n_channels>(std::forward<Args>(args)...);
            }
            return FactoryImpl<n_channels-1>::create(n, std::forward<Args>(args)...);
        }
    };

    // Base case specialization
    template<typename Dummy>
    struct FactoryImpl<0, Dummy> {
        template<typename... Args>
        static ISrc* create(int n_channels, Args&&... args) {
            throw std::invalid_argument("Invalid parameters");
        }
    };
}

class SrcFactory {
public:
    template<typename... Args>
    static ISrc* create(int n_channels, Args&&... args) {
        if (n_channels <= 0 || n_channels > 2) {
            throw std::invalid_argument("n_channels must be between 1 and 2");
        }
        return factory_detail::FactoryImpl<2>::create(n_channels, std::forward<Args>(args)...);
    }
};
