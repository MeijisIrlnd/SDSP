//
// Created by Syl on 29/08/2023.
//

#ifndef KALIDEPROTOTYPES2_SVFCOEFFS_H
#define KALIDEPROTOTYPES2_SVFCOEFFS_H
#include <SDSP/KMath.h>
#include <span.hpp>
#include <cmath>
#include <gcem.hpp>
namespace SDSP::Filters {
    enum class SVF_TYPE {
        HIGHPASS,
        BANDPASS,
        LOWPASS,
        NOTCH,
        PEAK,
        ALLPASS,
        BELL,
        LOWSHELF,
        HIGHSHELF,
    };
    enum class BELL_TYPE {
        Q,
        BW,
        SLOPE
    };
    template <SVF_TYPE SVFType>
    class SVFCoeffs {
    public:
        template <SVF_TYPE Type = SVFType, typename std::enable_if_t<Type == SVF_TYPE::HIGHPASS>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = -k;
            m_targetScalar[2] = -1.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::BANDPASS>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 0.0f;
            m_targetScalar[1] = 1.0f;
            m_targetScalar[2] = 0.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::LOWPASS>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 0.0f;
            m_targetScalar[1] = 0.0f;
            m_targetScalar[2] = 1.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::NOTCH>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = -k;
            m_targetScalar[2] = 0.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::PEAK>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = -k;
            m_targetScalar[2] = -2.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::ALLPASS>::value* = nullptr>
        void calculate(double sampleRate, float cutoff, float q) noexcept {
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / q;
            m_targetA[0] = 1 / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = -2.0f * k;
            m_targetScalar[2] = 0.0f;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::BELL>::type* = nullptr>
        void calculate(BELL_TYPE bellType, double sampleRate, float cutoff, float qualityParam, float dbGain) noexcept {
            const auto A = std::powf(10.0f, dbGain / 40.0f);
            float q;
            switch (bellType) {
                default: [[fallthrough]];
                case BELL_TYPE::Q: {
                    q = qualityParam;
                    break;
                }
                case BELL_TYPE::BW: {
                    static constexpr float ln2_over_2 = gcem::log<float>(2.0f) / 2.0f;
                    // static constexpr float ln2_over_2 = gcem::log(2) / 2.0f;
                    const auto omega = (static_cast<float>(m_pi) * 2.0f) * (cutoff / static_cast<float>(sampleRate));
                    const auto omega_over_sin_omega = omega / std::sinf(omega);
                    const auto inverse_q = 2.0f * std::sinh(ln2_over_2 * qualityParam * omega_over_sin_omega);
                    q = 1.0f / inverse_q;
                    break;
                }
                case BELL_TYPE::SLOPE: {
                    const auto inner1 = A + (1.0f / A);
                    const auto inner2 = (1.0f / qualityParam) - 1.0f;
                    const auto innerProduct = inner1 * inner2;
                    const auto inverseQ = std::sqrt(innerProduct + 2.0f);
                    q = 1.0f / inverseQ;
                    break;
                }
            }

            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate));
            const auto k = 1.0f / (q * A);
            m_targetA[0] = 1.0f / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = k * (A * A - 1);
            m_targetScalar[2] = 0.0f;
        }


        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::LOWSHELF>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q, float dbGain) noexcept {
            const auto A = std::powf(10.0f, dbGain / 40.0f);
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate)) / std::sqrt(A);
            const auto k = 1.0f / q;
            m_targetA[0] = 1.0f / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = 1.0f;
            m_targetScalar[1] = k * (A - 1);
            m_targetScalar[2] = A * A - 1;
        }

        template <SVF_TYPE Type = SVFType, typename std::enable_if<Type == SVF_TYPE::HIGHSHELF>::type* = nullptr>
        void calculate(double sampleRate, float cutoff, float q, float dbGain) noexcept {
            const auto A = std::powf(10.0f, dbGain / 40.0f);
            const auto g = std::tanf((static_cast<float>(m_pi) * cutoff) / static_cast<float>(sampleRate)) / std::sqrt(A);
            const auto k = 1.0f / q;
            m_targetA[0] = 1.0f / (1 + g * (g + k));
            m_targetA[1] = g * m_targetA[0];
            m_targetA[2] = g * m_targetA[1];
            m_targetScalar[0] = A * A;
            m_targetScalar[1] = k * (1 - A) * A;
            m_targetScalar[2] = 1 - A * A;
        }

        [[nodiscard]] tcb::span<float, 3> a() {
            return tcb::make_span(m_currentA);
        }

        [[nodiscard]] tcb::span<float, 3> scalar() {
            return tcb::make_span(m_currentScalar);
        }


        void interpolate() noexcept {
            for (size_t i = 0; i < 3; ++i) {
                m_currentA[i] = m_currentA[i] - 0.004f * (m_currentA[i] - m_targetA[i]);
                m_currentScalar[i] = m_currentScalar[i] - 0.004f * (m_currentScalar[i] - m_targetScalar[i]);
            }
        }

        template <SVF_TYPE T2>
        std::tuple<std::array<float, 3>, std::array<float, 3>> interpolateTo(const SVFCoeffs<T2>& dest, float proportion) {
            std::array<float, 3> resA{ 0.0f, 0.0f, 0.0f }, resScalar{ 0.0f, 0.0f, 0.0f };
            for (size_t i = 0; i < 3; ++i) {
                resA[i] = SDSP::KMath::Lerp<float>(m_targetA[i], dest.m_targetA[i], proportion);
                resScalar[i] = SDSP::KMath::Lerp<float>(m_targetScalar[i], m_targetScalar[i], proportion);
            }
            return { resA, resScalar };
        }

        [[nodiscard]] std::array<float, 3>& currentA() noexcept { return m_currentA; }
        [[nodiscard]] std::array<float, 3>& targetA() noexcept { return m_targetA; }
        [[nodiscard]] std::array<float, 3>& currentScalar() noexcept { return m_currentScalar; };
        [[nodiscard]] std::array<float, 3>& targetScalar() noexcept { return m_targetScalar; }

    private:
        std::array<float, 3> m_currentA{ 0.0f, 0.0f, 0.0f }, m_targetA{ 0.0f, 0.0f, 0.0f };
        std::array<float, 3> m_currentScalar{ 0.0f, 0.0f, 0.0f }, m_targetScalar{ 0.0f, 0.0f, 0.0f };
        constexpr static double m_pi = 3.14159265358979323846;
    };

} // namespace SDSP::Filters


#endif // KALIDEPROTOTYPES2_SVFCOEFFS_H
