#pragma once
#include <memory>

#include "converter.h"

namespace factory_detail
{
    // Helper template for generating the correct class based on parameters
    template<int n_channels, Quality quality, typename... Args>
    ISrc* createImpl(Args&&... args) {
        return new Src<n_channels, quality>(std::forward<Args>(args)...);
    }

    // Factory implementation that uses template recursion for channels
    template<int n_channels, typename = void>
    struct ChannelFactoryImpl {
        template<typename... Args>
        static ISrc* create(int n, Quality quality, Args&&... args)
        {
            if (n == n_channels) {
                return createQualityImpl<n_channels>(quality, std::forward<Args>(args)...);
            }
            return ChannelFactoryImpl<n_channels-1>::create(n, quality, std::forward<Args>(args)...);
        }

    private:
        template<int channels, typename... Args>
        static ISrc* createQualityImpl(Quality quality, Args&&... args) {
            switch (quality) {
                case Quality::LOW:
                    return createImpl<channels, Quality::LOW>(std::forward<Args>(args)...);
                case Quality::MEDIUM:
                    return createImpl<channels, Quality::MEDIUM>(std::forward<Args>(args)...);
                case Quality::HIGH:
                    return createImpl<channels, Quality::HIGH>(std::forward<Args>(args)...);
                default:
                    throw std::invalid_argument("Invalid quality parameter");
            }
        }
    };

    // Base case specialization for channels
    template<typename Dummy>
    struct ChannelFactoryImpl<0, Dummy> {
        template<typename... Args>
        static ISrc* create(int n_channels, Quality quality, Args&&... args) {
            throw std::invalid_argument("Invalid parameters");
        }
    };
}

class SrcFactory {
public:
    template<typename... Args>
    static ISrc* create(int n_channels, Quality quality, Args&&... args) {
        if (n_channels <= 0 || n_channels > 2) {
            throw std::invalid_argument("n_channels must be between 1 and 2");
        }
        return factory_detail::ChannelFactoryImpl<2>::create(n_channels, quality, std::forward<Args>(args)...);
    }
};