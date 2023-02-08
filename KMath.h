/*
  ==============================================================================

    KMath.h
    Created: 28 Dec 2021 4:04:30pm
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <inttypes.h>
#include <cmath>
#include "Macros.h"
namespace SDSP
{
    namespace KMath
    {
        SDSP_UNUSED static float fastLog2(float val) {
            //union { float val; int32_t x; } u = { val };
            //register float log2 = static_cast<float>(((u.x >> 23) & 255) - 128);
            //u.x &= ~(255 << 23);
            //u.x += 127 << 23;
            //log2 += ((-0.3358287811f) * u.val + 2.0f) * u.val - 0.65871759316667f;
            //return log2;
            return std::log2f(val);
        }

        static double fastPow(double a, double b) {
            //union {
            //    double d;
            //    int x[2];
            //} u = { a };
            //u.x[1] = static_cast<int>(b * (u.x[1] - 1072632447) + 1072632447);
            //u.x[0] = 0;
            //return u.d;
            return std::pow(a, b);
        }

        template<typename T, typename U>
        constexpr double dmod(T x, U mod) {
            return !mod ? x : x - mod * static_cast<long long>(x / mod);
        }

        template<typename T> 
        static T Lerp(T start, T end, float distance)
        {
            return static_cast<T>(start + (end - start) * distance);
        }

        template<typename T>
        static inline T log(T x, T base) { 
            
            auto logA = std::log2(x);
            auto logB = std::log2(base);
            return logA / logB;
        }
        static inline bool isPrime(int x) {
            for(int i = 3; i < std::sqrt(x); i += 2) {
                if(x % i == 0) {
                    return false;
                }
            }
            return true;
        }

        static inline float getNearestCoprime(float toCheck)
        {
            int above = std::ceil(toCheck);
            int below = std::floor(toCheck);
            if(above <= 2) {
                return 2;
            }
            if(below == 2) {
                return (toCheck - 2 < 0.5f) ? 2 : 3;
            }
            if(below % 2 == 0) {
                below -= 1;
            }
            if(above % 2 == 0) {
                above += 1;
            }

            double deltaBelow = std::numeric_limits<double>::max(), deltaAbove = std::numeric_limits<double>::max();;
            for(;; above += 2, below -= 2) {
                if(isPrime(below)){
                    deltaBelow = toCheck - below;
                }
                if(isPrime(above)) {
                    deltaAbove = above - toCheck;
                }
                if(deltaAbove != std::numeric_limits<double>::max() || deltaBelow != std::numeric_limits<double>::max()) {
                    break;
                }
            }
            return static_cast<int>(deltaAbove < deltaBelow ? above : below);
        }

        static inline std::vector<float> s_primes = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
            109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
            233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
            367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
            499, 503, 509, 521, 523, 541
        };

        static inline float nearestPrime(float p_i, float M_i) {
            auto m_i = std::roundf(std::logf(M_i) / std::logf(p_i));
            return std::powf(p_i, m_i);
        }

    }
}
