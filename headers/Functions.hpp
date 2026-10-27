#pragma once
#include <type_traits>
class Functions
{
    public:
        template<typename T>
        static auto normalCDF(T x) -> std::common_type_t <double , T>
        {
            using U = std::common_type_t<double, T>;
            return static_cast<U>(0.5) * (static_cast<U>(1) + std::erf(static_cast<U>(x) / std::sqrt(static_cast<U>(2.0))));
        }
};