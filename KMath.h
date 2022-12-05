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
        static T Lerp(T start, T end, double distance)
        {
            return static_cast<T>(start + (end - start) * distance);
        }
    }
}
