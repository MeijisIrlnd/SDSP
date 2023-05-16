//
// Created by Syl on 04/02/2023.
// Potentially not that useful, but nice drop ins to try and avoid using std::mem stuff directly
//
#pragma once
#include "../Macros.h"
#include <array>
#include <concepts>
namespace SDSP::Helpers
{
    template<class T>
    struct is_std_array : std::is_array<T> {};

    template<class T, size_t N>
    struct is_std_array<typename std::array<T, N> > : std::true_type { };

#if __cplusplus >= 202002L
    template<typename T>
    concept NumericArray =
        (SDSP::Helpers::is_std_array<T>::value && std::is_arithmetic<typename T::value_type>::value) || std::is_array<T>::value;
#endif
#if __cplusplus >= 202002L
    template<typename T> requires NumericArray<T>
#else
    template<typename T>
#endif
    SDSP_UNUSED static inline void zero_array(T& data) {
        std::memset(data.data(), 0.0f, sizeof(typename T::value_type) * data.size());
    }
#if __cplusplus >= 202002L
    template<typename T> requires NumericArray<T>
#else
    template<typename T>
#endif
    SDSP_UNUSED static inline void copy_array(T& dest, T& src) {
        std::memcpy(dest.data(), src.data(), sizeof(typename T::value_type) * dest.size());
    }
#if __cplusplus >= 202002L
    template<typename T> requires NumericArray<T>
#else
    template<typename T>
#endif
    SDSP_UNUSED static inline void fill_array(T& dest, typename T::value_type value) {
        std::fill(dest.begin(), dest.end(), value);
    }
}