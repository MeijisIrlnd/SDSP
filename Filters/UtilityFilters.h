//
// Created by Syl on 05/03/2023.
//
#pragma once
#include "DSPBiquad.h"
#include "RBJCoefficients.h"
#include "SmoothedFilterCoefficients.h"
namespace SDSP::Filters
{
    class [[maybe_unused]] SinglePoleLowpass
    {
    public:

        void setCoeff(const float newCoeff) noexcept {
            m_coeff = newCoeff;
        }

        float processSample(float in) {
            in *= 1 - m_coeff;
            auto current = in + m_x1;
            m_x1 = current * m_coeff;
            return current;
        }


    private:
        float m_x1{ 0.0f };
        float m_coeff{ 0.9995f };
    };

    class TiltShift {
    public:
        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate) {
            m_sampleRate = sampleRate;
            SDSP::RBJ::lowShelf(m_lowShelfCoeffs.target(0), sampleRate, m_xOverFreq, m_lowShelfGainDb, m_slope);
            std::memcpy(m_lowShelfCoeffs.current(0), m_lowShelfCoeffs.target(0), sizeof(double) * 6);
            SDSP::RBJ::highShelf(m_highShelfCoeffs.target(0), sampleRate, m_xOverFreq, -1.0f * m_lowShelfGainDb, m_slope);
            std::memcpy(m_highShelfCoeffs.current(0), m_highShelfCoeffs.target(0), sizeof(double) * 6);
            m_lowShelf.setCoefficients(m_lowShelfCoeffs.target(0));
            m_highShelf.setCoefficients(m_highShelfCoeffs.target(0));
        }

        float processSample(float x) {
            if(m_samplesUntilUpdate == 0) {
                SDSP::RBJ::lowShelf(m_lowShelfCoeffs.target(0), m_sampleRate, m_xOverFreq, m_lowShelfGainDb, m_slope);
                SDSP::RBJ::highShelf(m_highShelfCoeffs.target(0), m_sampleRate, m_xOverFreq, -1.0f * m_lowShelfGainDb, m_slope);
                m_samplesUntilUpdate = m_updateRate;
            }
            interpolate();
            --m_samplesUntilUpdate;
            x = m_lowShelf.processSample(x);
            x = m_highShelf.processSample(x);
            return x;
        }

        SDSP_INLINE void setCrossoverFreq(float newCrossover) {
            m_xOverFreq = newCrossover;
        }

        SDSP_INLINE void setLowShelfGain(float newGain) {
            m_lowShelfGainDb = newGain;
        }

        SDSP_INLINE void setSlope(float newSlope) {
            m_slope = newSlope;
        }
    private:
        void interpolate() {
            m_lowShelfCoeffs.interpolate();
            m_highShelfCoeffs.interpolate();
            m_lowShelf.setCoefficients(m_lowShelfCoeffs.current(0));
            m_highShelf.setCoefficients(m_highShelfCoeffs.current(0));
        }

        const int m_updateRate{ 100 };
        int m_samplesUntilUpdate{ 0 };
        double m_sampleRate{};
        BiquadCascade<1> m_lowShelf, m_highShelf;
        SmoothedFilterCoefficients<1> m_lowShelfCoeffs, m_highShelfCoeffs;
        float m_xOverFreq{ 1000.0f };
        float m_lowShelfGainDb{ 0.0f };
        float m_slope{ 0.5f };
    };
}
