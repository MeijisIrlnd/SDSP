//
// Created by Syl on 04/02/2023.
// Potentially not that useful, but nice drop ins to try and avoid using std::mem stuff directly
//
#pragma once
#include "../Macros.h"
#include <array>
#include <concepts>
#include <juce_core/juce_core.h>
namespace SDSP::Helpers {
    template <class T>
    struct is_std_array : std::is_array<T> {};

    template <class T, size_t N>
    struct is_std_array<typename std::array<T, N>> : std::true_type {};

#if __cplusplus >= 202002L
    template <typename T>
    concept NumericArray =
        (SDSP::Helpers::is_std_array<T>::value && std::is_arithmetic<typename T::value_type>::value) || std::is_array<T>::value;
#endif
#if __cplusplus >= 202002L
    template <typename T>
    requires NumericArray<T>
#else
    template <typename T>
#endif
        SDSP_UNUSED static inline void zero_array(T& data) {
        std::memset(data.data(), 0.0f, sizeof(typename T::value_type) * data.size());
    }
#if __cplusplus >= 202002L
    template <typename T>
    requires NumericArray<T>
#else
    template <typename T>
#endif
        SDSP_UNUSED static inline void copy_array(T& dest, T& src) {
        std::memcpy(dest.data(), src.data(), sizeof(typename T::value_type) * dest.size());
    }
#if __cplusplus >= 202002L
    template <typename T>
    requires NumericArray<T>
#else
    template <typename T>
#endif
        SDSP_UNUSED static inline void fill_array(T& dest, typename T::value_type value) {
        std::fill(dest.begin(), dest.end(), value);
    }

    inline void protectYourEars(float* buffer, int sampleCount) {
        if (buffer == nullptr) {
            return;
        }
        bool firstWarning = true;
        for (int i = 0; i < sampleCount; ++i) {
            float x = buffer[i];
            bool silence = false;
            if (std::isnan(x)) {
                DBG("!!! WARNING: nan detected in audio buffer, silencing !!!");
                silence = true;
            } else if (std::isinf(x)) {
                DBG("!!! WARNING: inf detected in audio buffer, silencing !!!");
                silence = true;
            } else if (x < -2.0f || x > 2.0f) { // screaming feedback
                DBG("!!! WARNING: sample out of range, silencing !!!");
                silence = true;
            } else if (x < -1.0f) {
                if (firstWarning) {
                    DBG("!!! WARNING: sample out of range, clamping !!!");
                    firstWarning = false;
                }
                buffer[i] = -1.0f;
            } else if (x > 1.0f) {
                if (firstWarning) {
                    DBG("!!! WARNING: sample out of range, clamping !!!");
                    firstWarning = false;
                }
                buffer[i] = 1.0f;
            }
            if (silence) {
                memset(buffer, 0, sampleCount * sizeof(float));
                return;
            }
        }
    }
} // namespace SDSP::Helpers