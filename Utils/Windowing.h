#pragma once
#include <juce_core/juce_core.h>
#include <juce_dsp/juce_dsp.h>
#if PERFETTO
    #include <melatonin_perfetto/melatonin_perfetto.h>
#endif
namespace SDSP::Utils::Windowing {

    template<typename ReturnType, typename IndexingType>
    [[maybe_unused]] [[nodiscard]] static inline ReturnType hann(IndexingType i, IndexingType n) {
#if PERFETTO
        TRACE_DSP();
#endif
        auto pi_n = juce::MathConstants<ReturnType>::pi * static_cast<ReturnType>(i);
        auto pi_n_over_N = pi_n / static_cast<ReturnType>(n);
        auto sin = juce::dsp::FastMathApproximations::sin<ReturnType>(pi_n_over_N);
        return sin * sin;
    }

    template<typename ReturnType, typename IndexingType>
    [[maybe_unused]] [[nodiscard]] static inline ReturnType tukey(IndexingType i, IndexingType n, ReturnType alpha) {
        auto lobeWidth = alpha * static_cast<ReturnType>(n);
        lobeWidth /= static_cast<ReturnType>(2);
        if(static_cast<ReturnType>(i) < lobeWidth) {
            auto angle = juce::MathConstants<ReturnType>::twoPi * static_cast<ReturnType>(i);
            angle /= alpha * static_cast<ReturnType>(n);
            auto cos = juce::dsp::FastMathApproximations::cos<ReturnType>(angle);
            auto inv = 1 - cos;
            inv *= static_cast<ReturnType>(0.5);
            return inv;
        }
        else if(lobeWidth <= static_cast<ReturnType>(i) && i <= static_cast<ReturnType>(n) - lobeWidth) {
            return 1;
        }
        else {
            auto invIndex = n - i;
            return tukey(invIndex, n, alpha);
        }

    }
}