//
// Created by Syl on 14/03/2023.
//
#pragma once
#include <juce_dsp/juce_dsp.h>
#include "../Macros.h"
#include "../KMath.h"
namespace SDSP
{
    class BlendableLFO
    {
    public:
        BlendableLFO(const std::function<float(float)>& generator1, const std::function<float(float)>& generator2) :
        m_generator1(generator1), m_generator2(generator2)
        {

        }

        SDSP_INLINE void setRate(float newRate) noexcept {
            m_rate = newRate;
            m_osc1.setFrequency(m_rate);
            m_osc2.setFrequency(m_rate);
        }

        SDSP_INLINE void setBlendPosition(float newBlendPos) noexcept {
            m_blendPos = newBlendPos;
        }

        void prepareToPlay(int samplesPerBlockExpected, double sampleRate) {
            juce::dsp::ProcessSpec spec{sampleRate, static_cast<juce::uint32>(samplesPerBlockExpected), 1};
            m_osc1.prepare(spec);
            m_osc1.initialise(m_generator1);

            m_osc2.prepare(spec);
            m_osc2.initialise(m_generator2);

        }

        [[nodiscard]] float getNextSample() noexcept {
            // 0 - 0.5 => 1 to 2
            // 0.5 - 1 => 2 to 3
            auto out1 = m_osc1.processSample(0.0f);
            auto out2 = m_osc2.processSample(0.0f);
            // ((max - min) * T) + min;
            auto blended = KMath::Lerp<float>(out1, out2, m_blendPos);
            return blended;
        }

    private:
        float m_blendPos{ 0.0f };
        float m_rate{ 0.01f };
        std::function<float(float)> m_generator1, m_generator2;
        juce::dsp::Oscillator<float> m_osc1, m_osc2;
    };
}