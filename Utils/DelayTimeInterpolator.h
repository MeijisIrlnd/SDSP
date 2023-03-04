//
// Created by Syl on 04/02/2023.
//
#pragma once
#include <array>
#include "../Helpers/Helpers.h"
namespace SDSP
{
    template<typename ArrayType> requires Helpers::NumericArray<ArrayType>
    struct SDSP_UNUSED DelayTimeInterpolator
    {
        DelayTimeInterpolator() {
            Helpers::zero_array(target);
            Helpers::zero_array(current);
            Helpers::zero_array(targetModulated);
            Helpers::zero_array(currentModulated);
        }

        void interpolate() {
            for(size_t i = 0; i < current.size(); i++) {
                current[i] = current[i] - 0.004 * (current[i] - target[i]);
                currentModulated[i] = currentModulated[i] - 0.004 * (currentModulated[i] - targetModulated[i]);
            }
        }

        void setCurrentAndTargetValue(size_t index, typename ArrayType::value_type value) noexcept {
            target[index] = current[index] = targetModulated[index] = current[index] = value;
        }

        ArrayType target, current, targetModulated, currentModulated;
    };
}