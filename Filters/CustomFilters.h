//
// Created by Syl on 07/02/2023.
//
#pragma once
#include <tuple>
#include "../Helpers/Helpers.h"
#include "../Macros.h"

namespace SDSP::Filters
{
    SDSP_UNUSED static inline void fdnDecayFilter(std::array<double, 2>& targetB, std::array<double, 2>& targetA, double sampleRate, int delayLineLength, double dcRT60, double halfSRRT60)
    {
        // g_i = 10^(-3MiT/t60(0))
        auto samplingPeriod = 1 / sampleRate;
        auto t0 = samplingPeriod / dcRT60;
        auto exponent = -3 * delayLineLength * t0;
        auto g_i = std::pow(10, exponent);
        //a = rt60(w / t) / rt60(0)
        auto alpha = halfSRRT60 / dcRT60;
        // p_i = (ln(10) / 4) * log_10(g_i)(1 - 1 / a^2)
        auto alphaExpression = 1 - (1 / std::pow(alpha, 2));
        auto ln10 = std::log(10) / 4.0;
        auto log10gi = std::log10(g_i);
        auto p_i = ln10 * log10gi * alphaExpression;

        targetB[0] = g_i * (1 - p_i);
        targetA[0] = 1;
        targetA[1] = -p_i;
    }

    SDSP_UNUSED static inline void fdnTonalCorrectionFilter(std::array<double, 2>& targetB, std::array<double, 2>& targetA, double dcRT60, double halfSRRT60)
    {
        auto alpha = halfSRRT60 / dcRT60;
        auto x = (1 - alpha) / (1 + alpha);
        targetB[0] = 1 / (1 - x);
        targetB[1] = -x / (1 - x);
        targetA[0] = 1;
    }

    SDSP_UNUSED static inline std::tuple<double, double> fdnDecayFilter(double sampleRate, int delayLineLength, double dcRT60, double halfSRRT60)
    {
        // g_i = 10^(-3MiT/t60(0))
        auto samplingPeriod = 1 / sampleRate;
        auto t0 = samplingPeriod / dcRT60;
        auto exponent = -3 * delayLineLength * t0;
        auto g_i = std::pow(10, exponent);
        //a = rt60(w / t) / rt60(0)
        auto alpha = halfSRRT60 / dcRT60;
        // p_i = (ln(10) / 4) * log_10(g_i)(1 - 1 / a^2)
        auto alphaExpression = 1 - (1 / std::pow(alpha, 2));
        auto ln10 = std::log(10) / 4.0;
        auto log10gi = std::log10(g_i);
        auto p_i = ln10 * log10gi * alphaExpression;
        return {(g_i * (1 - p_i)), -p_i};
    }

    struct SDSP_UNUSED SinglePole
    {
        SDSP_INLINE float processSample(float in) noexcept {
            // y(n) = b0x(n) - a1x(n-1)
            coeffs.interpolate();
            float x = static_cast<float>(coeffs.currentBCoeffs[0] * in) + static_cast<float>(coeffs.currentBCoeffs[1] * m_prev) - static_cast<float>(coeffs.currentACoeffs[1] * m_prev);
            m_prev = x;
            return x;
        }

        struct Coeffs {
            std::array<double, 2> currentBCoeffs = {0.0f, 0.0f}, targetBCoeffs = {0.0f, 0.0f};
            std::array<double, 2> currentACoeffs = {1.0f, 0.0f}, targetACoeffs = {1.0f, 0.0f};

            SDSP_UNUSED SDSP_INLINE void setCurrentAndTarget(const std::array<double, 2>& bCoeffs, const std::array<double, 2>& aCoeffs) noexcept {
                juce::FloatVectorOperations::copy(targetBCoeffs.data(), bCoeffs.data(), 2);
                juce::FloatVectorOperations::copy(currentBCoeffs.data(), bCoeffs.data(), 2);
                juce::FloatVectorOperations::copy(targetACoeffs.data(), aCoeffs.data(), 2);
                juce::FloatVectorOperations::copy(currentACoeffs.data(), aCoeffs.data(), 2);
            }

            SDSP_UNUSED SDSP_INLINE void setTarget(const std::array<double, 2>& bCoeffs, const std::array<double, 2>& aCoeffs) noexcept {
                juce::FloatVectorOperations::copy(targetBCoeffs.data(), bCoeffs.data(), 2);
                juce::FloatVectorOperations::copy(targetACoeffs.data(), aCoeffs.data(), 2);
            }

            SDSP_UNUSED SDSP_INLINE void interpolate() noexcept {
                for(size_t i = 0; i < 2; ++i) {
                    currentBCoeffs[i] = currentBCoeffs[i] - 0.004 * (currentBCoeffs[i] - targetBCoeffs[i]);
                    currentACoeffs[i] = currentACoeffs[i] - 0.004 * (currentACoeffs[i] - targetACoeffs[i]);
                }
            }
        } coeffs;
    private:
        float m_prev{};
    };
}