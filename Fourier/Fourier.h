#pragma once 
#include "STFT.h"
namespace SDSP::Fourier
{
    template <typename T>
    static inline T mag(T real, T imag) {
        return std::sqrtf(real * real + imag * imag);
    }

    template <typename T>
    static inline T phase(T real, T imag) {
        return std::atan2<T>(imag, real);
    }

    // Principal argument - Unwrap a phase argument to between [-PI, PI]
    static inline float principalArgument(float arg)
    {
        return std::fmod(arg + juce::MathConstants<float>::pi,
            -juce::MathConstants<float>::twoPi) + juce::MathConstants<float>::pi;
    }
}