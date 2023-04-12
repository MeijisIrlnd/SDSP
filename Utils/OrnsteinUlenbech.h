//
// Created by Syl on 12/04/2023.
//

#pragma once
#include "../Oscillators/SDSPOscillator.h"
#include "../Filters/RBJCoefficients.h"
#include "../Filters/DSPBiquad.h"
#include "../Filters/SmoothedFilterCoefficients.h"
namespace SDSP {
    class OrnsteinUlenbech {
    public:
        OrnsteinUlenbech() : m_noiseGen(false, Oscillators::SHAPE::WHITE_NOISE) {

        }

        OrnsteinUlenbech(const OrnsteinUlenbech & /*other*/) : m_noiseGen(false, Oscillators::SHAPE::WHITE_NOISE) {}

        OrnsteinUlenbech(OrnsteinUlenbech && /*other*/) noexcept: m_noiseGen(false, Oscillators::SHAPE::WHITE_NOISE) {}

        void prepare(int samplesPerBlockExpected, double sampleRate) {
            m_noiseGen.prepareToPlay(samplesPerBlockExpected, sampleRate);
            m_sqrtDelta = 1.0f / std::sqrtf(static_cast<float>(sampleRate));
            m_time = 1.0f / static_cast<float>(sampleRate);
            RBJ::lowpass(m_coeffs.target(0), sampleRate, 10.0, 0.5);
            std::memcpy(m_coeffs.current(0), m_coeffs.target(0), sizeof(double) * 6);
            m_lowpass.setCoefficients(m_coeffs.target(0));
        }

        float processSample() {
            // Gain found here: https://github.com/jatinchowdhury18/AnalogTapeModel/blob/master/Plugin/Source/Processors/Timing_Effects/OHProcess.h
            float noise = (m_noiseGen.processSample() * (1.0f / 2.33f));
            m_prevSample += m_sqrtDelta * noise * m_amt;
            m_prevSample += m_damping * (m_mean - m_prevSample) * m_time;
            return m_lowpass.processSample(m_prevSample);
        }

        void setAmt(float newAmt) {
            m_amt = std::powf(newAmt, 1.25f);
            m_damping = m_amt * 20.0f + 1.0f;
            m_mean = newAmt;
        }

    private:
        SDSP::Oscillators::SDSPOscillator m_noiseGen;
        float m_sqrtDelta{}, m_time{};
        float m_amt{0.0f}, m_mean{0.0f}, m_damping{0.0f};
        BiquadCascade<1> m_lowpass;
        SmoothedFilterCoefficients<1> m_coeffs{};

        float m_prevSample{1.0f};
    };
}