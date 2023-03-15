//
// Created by Syl on 14/03/2023.
//
#pragma once
#include <juce_dsp/juce_dsp.h>
#include "../Macros.h"
#include "../KMath.h"
#include <SDSP/Oscillators/SDSPOscillator.h>
namespace SDSP
{
    class BlendableLFO
    {
    public:
        BlendableLFO() = default;

        BlendableLFO(SDSP::Oscillators::SHAPE shape1, SDSP::Oscillators::SHAPE shape2) : m_osc1(true), m_osc2(true) {
            m_osc1.setShape(shape1);
            m_osc2.setShape(shape2);
        }

        [[maybe_unused]] void setShapes(SDSP::Oscillators::SHAPE shape1, SDSP::Oscillators::SHAPE shape2) {
            m_osc1.setShape(shape1);
            m_osc2.setShape(shape2);
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
            m_osc1.prepareToPlay(samplesPerBlockExpected, sampleRate);
            m_osc2.prepareToPlay(samplesPerBlockExpected, sampleRate);
            m_osc1.setFrequency(m_rate);
            m_osc2.setFrequency(m_rate);
            m_osc1.retrigger();
            m_osc2.retrigger();
        }

        [[nodiscard]] float getNextSample() noexcept {
            // 0 - 0.5 => 1 to 2
            // 0.5 - 1 => 2 to 3
            auto out1 = m_osc1.processSample();
            auto out2 = m_osc2.processSample();
            // ((max - min) * T) + min;
            auto blended = KMath::Lerp<float>(out1, out2, m_blendPos);
            return blended;

        }

    private:
        double m_sampleRate{ 44100 };
        float m_blendPos{ 0.0f };
        float m_rate{ 0.01f };
        std::function<float(float)> m_generator1, m_generator2;
        Oscillators::SDSPOscillator m_osc1, m_osc2;
    };
}