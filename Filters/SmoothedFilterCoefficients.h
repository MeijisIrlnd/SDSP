/*
  ==============================================================================

    SmoothedFilterCoefficients.h
    Created: 3 Dec 2022 5:26:28pm
    Author:  Syl

  ==============================================================================
*/

#pragma once

namespace SDSP
{
    template<int NSTAGES>
    struct SmoothedFilterCoefficients
    {
        std::array<std::array<double, 6>, NSTAGES> currents;
        std::array<std::array<double, 6>, NSTAGES> targets;
        double* current(int stage) {
            return currents[stage].data();
        }

        double* target(int stage) {
            return targets[stage].data();
        }

        void interpolate()
        {
            for (auto stage = 0; stage < NSTAGES; stage++) {
                for (auto i = 0; i < 6; i++) {
                    currents[stage][i] = currents[stage][i] - 0.004 * (currents[stage][i] - targets[stage][i]);
                }
            }
        }
    };
}