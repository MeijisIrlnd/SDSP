//
// Created by Syl on 12/04/2023.
//
#pragma once
#include <array>

namespace SDSP::Filters {
//https://www.musicdsp.org/en/latest/Filters/85-1st-and-2nd-order-pink-noise-filters.html
    constexpr static inline std::array<double, 6> pinking() {
        return std::array<double, 6>{0.04957526213389, -0.06305581334498, 0.01483220320740, 1.00000000000000, -1.80116083982126, 0.80257737639225};
    }
}

