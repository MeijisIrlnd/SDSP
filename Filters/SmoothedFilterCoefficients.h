/*
  ==============================================================================

    SmoothedFilterCoefficients.h
    Created: 3 Dec 2022 5:26:28pm
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <array>
namespace SDSP
{
    template<int NSTAGES>
    struct [[maybe_unused]] SmoothedFilterCoefficients
    {
    private:
        std::array<double, 6> m_temp;
    public:

        std::array<std::array<double, 6>, NSTAGES> currents;
        std::array<std::array<double, 6>, NSTAGES> targets;
        [[nodiscard]] double* current(int stage) {
            return currents[static_cast<size_t>(stage)].data();
        }

        [[nodiscard]] double* target(int stage) {
            return targets[static_cast<size_t>(stage)].data();
        }

        void interpolate()
        {
            for(size_t stage = 0; stage < NSTAGES; stage++) {
                for (size_t i = 0; i < 6; i++) {
                    currents[stage][i] = currents[stage][i] - 0.004 * (currents[stage][i] - targets[stage][i]);
                }
            }
        }
    };
}