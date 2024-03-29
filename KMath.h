/*
  ==============================================================================

    KMath.h
    Created: 28 Dec 2021 4:04:30pm
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <vector>
#include <cinttypes>
#include <cmath>
#if PERFETTO
/*#include <melatonin_perfetto/melatonin_perfetto.h>*/
#endif

#if defined(FAST_MATH)
// Keep this around for when we inevitably want to increase perf...
#include <fastmaths/pow.hpp>
#endif
#include "Macros.h"
namespace SDSP::KMath {
    template <typename T>
    SDSP_UNUSED static T Lerp(T start, T end, float distance) {
        distance = distance > 1 ? 1 : distance;
        distance = distance < 0 ? 0 : distance;
        return static_cast<T>(start + (end - start) * distance);
    }

    template <typename T>
    SDSP_UNUSED static inline T log(T x, T base) {
        if (x == 0) throw std::exception();
        if (base == 0) return 0;
        auto logA = std::log2(x);
        auto logB = std::log2(base);
        return logA / logB;
    }

    SDSP_UNUSED static inline bool isPrime(int x) {
        if (x < 2) return false;
        for (auto i = 2; i <= x / 2; ++i) {
            if (x % i == 0) return false;
        }
        return true;
    }

    // Is this actually just nearest prime?
    SDSP_UNUSED static inline float getNearestCoprime(float toCheck) {
        const auto dToCheck{ static_cast<double>(toCheck)};
        auto above = static_cast<int>(std::ceil(toCheck));
        auto below = static_cast<int>(std::floor(toCheck));
        if (above <= 2) {
            return 2;
        }
        if (below == 2) {
            return (toCheck - 2.0f < 0.5f) ? 2.0f : 3.0f;
        }
        if (below % 2 == 0) {
            below -= 1;
        }
        if (above % 2 == 0) {
            above += 1;
        }

        double deltaBelow = std::numeric_limits<double>::max(), deltaAbove = std::numeric_limits<double>::max();
        ;
        for (;; above += 2, below -= 2) {
            if (isPrime(below)) {
                deltaBelow = dToCheck - static_cast<double>(below);
            }
            if (isPrime(above)) {
                deltaAbove = static_cast<double>(above) - dToCheck;
            }
            const auto diffAbv{ std::numeric_limits<double>::max() - std::abs(deltaAbove) };
            const auto diffBelow{ std::numeric_limits<double>::max() - std::abs(deltaBelow)};
            const auto epsilon{ std::numeric_limits<double>::min() };
            if(diffAbv > epsilon || diffBelow > epsilon) {
                break;
            }
/*
            if (deltaAbove != std::numeric_limits<double>::max() || deltaBelow != std::numeric_limits<double>::max()) {
                break;
            }
*/
        }
        return static_cast<float>(deltaAbove < deltaBelow ? above : below);
    }

    static inline std::vector<float> s_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541
    };

    [[maybe_unused]] [[nodiscard]] static inline float nearestPrime(float source) {
        if (source < 2) return 2;
        auto res = source;
        while (!isPrime(static_cast<int>(res))) {
            --res;
        }
        return res;
    }

    // TODO: UNIT TEST ME
    [[maybe_unused]] [[nodiscard]] static inline double polyBlep(double t, double phaseIncrement) noexcept {
        // this gets you out of radians, we're already out of radians, with 0 to 1 phase
        // double dt = phaseIncrement / juce::MathConstants<float>::twoPi;
        double dt = phaseIncrement;
        if (t < dt) {
            t /= dt;
            return t + t - t * t - 1.0;
        } else if (t > 1.0 - dt) {
            t = (t - 1.0) / dt;
            return t * t + t + t + 1.0;
        }
        return 0;
    }

    [[maybe_unused]] [[nodiscard]] static bool usingIPP() {
#if defined(HAS_IPP)
        return true;
#else
        return false;
#endif
    }


    template <typename Type>
    [[maybe_unused]] [[nodiscard]] static Type decibelsToGain(Type decibels, Type minusInfinityDb = Type(-100)) {
        return decibels > minusInfinityDb ? pow<Type>(static_cast<Type>(10.0), decibels * static_cast<Type>(0.05)) : Type();
    }
} // namespace SDSP::KMath
