/*
  ==============================================================================

    Biquad.h
    Created: 1 Jan 2022 5:47:04pm
    Author:  Syl

  ==============================================================================
*/

#pragma once

#include <cinttypes>
#include <vector>
#include <juce_core/juce_core.h>
#include <juce_dsp/juce_dsp.h>
#include "../Macros.h"
#include "RBJCoefficients.h"
namespace SDSP
{
    template<int N>
    constexpr static inline std::array<double, N / 2> getButterworthQs()
    {
        std::array<double, N / 2> res;
        constexpr auto poleSpacing = juce::MathConstants<double>::pi / static_cast<double>(N);
        if constexpr(N % 2 == 0) {
            // First pole is poleSpacing /  2 radians from the real axis
            for (auto i = 0; i < N / 2; i++) {
                res[i] = 1 / (2 * juce::dsp::FastMathApproximations::cos<double>((poleSpacing / 2.0) + (i * poleSpacing)));
            }
        }
        else {
            for (auto i = 0; i < N / 2; i++) {
                // spacing + i * spacing right? 
                res[i] = 1 / (2 * juce::dsp::FastMathApproximations::cos<double>(poleSpacing + (i * poleSpacing)));
            }
        }
        return res;
    }

    class BiQuadDelay
    {
    public:
        BiQuadDelay() = default;
        BiQuadDelay(const BiQuadDelay& other) = default;

        float mX_z1;
        float mX_z2;
        float mY_z1;
        float mY_z2;

        void reset() {
            mX_z1 = mX_z2 = mY_z1 = mY_z2 = 0.0;
        }

    };

    class Biquad
    {
    public:
        Biquad(const uint32_t numSections, SDSP_UNUSED const uint32_t maxFramesPerSlice) : m_numSections(numSections)
        {
            m_biquadDelays.resize(m_numSections);
            m_coefficients = std::vector<float>(m_numSections * 6, 0.0);
            for (uint32_t sec = 0; sec < numSections; sec++) {
                m_coefficients.at(0 + 6 * sec) = 1.0;
            }
            reset();
        }

        Biquad(const Biquad& other) : m_numSections(other.m_numSections), m_biquadDelays(other.m_biquadDelays), m_coefficients(other.m_coefficients)
        {
            reset();
        }

        void initialise(const double* coeffs) {
            for (uint32_t i = 0; i < m_numSections * 6; i++) {
                m_coefficients[i] = static_cast<float>(coeffs[i]);
            }
        }

        template<int NumFrames>
        void processSample(const float* input, float* output) {
            //    //v(n) = gx(n)
            // y(n) = B1v(n - 1) + B2v(n - 2) - a1(
            // If we write dif eq as y = a0*x + a1*xz_1 + a2*xz_2 - b0*yz_1 - b1*yz_2, 
            // 
            //then the passed coeffs arg should be a concatenated groups of 5 coeffs containing 
            //{a0, a1, a2, b0, b1}, with the number of these groups determining the number of sections
            for (uint32_t frame = 0; frame < NumFrames; frame++) {
                float sectionInput = input[frame];
                for (uint32_t sec = 0; sec < m_numSections; sec++) {
                    uint32_t offset = sec * 5;
                    BiQuadDelay* del = &m_biquadDelays[sec];
                    float y = 1 / m_coefficients[offset + 3] * ( // b0
                        m_coefficients[offset] * sectionInput + // a0
                        m_coefficients[offset + 1] * del->mX_z1 + // a1
                        m_coefficients[offset + 2] * del->mX_z2 - // a2
                        m_coefficients[offset + 4] * del->mY_z1 - // b1
                        m_coefficients[offset + 5] * del->mY_z2  // b2
                        );
                    del->mX_z2 = del->mX_z1;
                    del->mX_z1 = sectionInput;
                    del->mY_z2 = del->mY_z1;
                    del->mY_z1 = y;

                    sectionInput = y;
                }
                output[frame] = sectionInput;
            }
        }
        void reset()
        {
            for (uint32_t i = 0; i < m_numSections; i++) {
                m_biquadDelays.at(i).reset();
            }
        }
        SDSP_UNUSED SDSP_NODISCARD uint32_t getNumSections() const { return m_numSections; }
    private:
        uint32_t m_numSections;
        std::vector<BiQuadDelay> m_biquadDelays;
        std::vector<float> m_coefficients;
    private:

    };
    
    template<int N>
    class BiquadCascade
    {
    public: 
        BiquadCascade() {
            for (auto i = 0; i < N; i++) {
                // m_biquads.add(new Biquad(1, 1));
                m_biquads.emplace_back(1, 1);
            }
        }

        float processSample(float in) {
            for (auto i = 0; i < N; i++) {
              m_biquads[i].template processSample<1>(&in, &in);
            }
            return in;
        }

        void setCoefficients(double* coeffs, int stage) {
            if (stage < N) {
                m_biquads[stage].initialise(coeffs);
            }
        }

        void setCoefficients(double* coeffs) {
            for (auto i = 0; i < N; i++) {
                m_biquads[i].initialise(coeffs);
            }
        }

        void reset() {
            for (auto& b : m_biquads) { b.reset(); }
        }

    private: 
        std::vector<Biquad> m_biquads;
    };
}
