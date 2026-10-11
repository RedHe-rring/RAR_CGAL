#pragma once

#include <array>

namespace rar {

using CSFMultipliers = std::array<double, 5>;

inline constexpr CSFMultipliers kDefaultCSFMultipliers{
    1.8, 1.4, 1.0, 0.8, 0.6};

} // namespace rar
