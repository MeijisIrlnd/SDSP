//
// Created by Syl on 28/08/2023.
//

#ifndef KALIDEPROTOTYPES2_CYTOMICTPT_H
#define KALIDEPROTOTYPES2_CYTOMICTPT_H
#include <span.hpp>
#include <array>
#include <atomic>
namespace SDSP::Filters {
    class CytomicTPT {
    public:
        [[maybe_unused]] [[nodiscard]] float processSample(float x) noexcept {
            const auto v3 = x - m_ic2Eq;
            const auto v1 = m_coeffs.a(0) * m_ic1Eq + m_coeffs.a(1) * v3;
            const auto v2 = m_ic2Eq + m_coeffs.a(1) * m_ic1Eq + m_coeffs.a(2) * v3;
            m_ic1Eq = 2.0f * v1 - m_ic1Eq;
            m_ic2Eq = 2.0f * v2 - m_ic2Eq;
            const auto output = m_coeffs.scalar(0) * x + m_coeffs.scalar(1) * v1 + m_coeffs.scalar(2) * v2;
            return output;
        }

        void setCoeffs(tcb::span<float, 3> aCoeffs, tcb::span<float, 3> scalars) noexcept {
            m_coeffs.setA(aCoeffs);
            m_coeffs.setScalars(scalars);
        }

    private:
        struct InternalCoeffs {
            std::array<float, 3> aCoeffs{ 0.0f, 0.0f, 0.0f };
            std::array<float, 3> scalars{ 0.0f, 0.0f, 0.0f };
            [[nodiscard]] float a(size_t index) noexcept { return aCoeffs[index]; }
            [[nodiscard]] float scalar(size_t index) noexcept { return scalars[index]; }
            void setA(tcb::span<float, 3> toSet) noexcept {
                std::memcpy(aCoeffs.data(), toSet.data(), 3 * sizeof(float));
            }
            void setScalars(tcb::span<float, 3> toSet) noexcept {
                std::memcpy(scalars.data(), toSet.data(), 3 * sizeof(float));
            }
        } m_coeffs;
        float m_ic1Eq{ 0.0f }, m_ic2Eq{ 0.0f };
    };
}


#endif //KALIDEPROTOTYPES2_CYTOMICTPT_H
