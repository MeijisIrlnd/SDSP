//
// Created by Syl on 15/03/2023.
//
#pragma once

#include "../Filters/DSPBiquad.h"
#include "../Filters/FilterHelpers.h"
#include "../Macros.h"

#include <juce_core/juce_core.h>

#include <cmath>

namespace SDSP::Oscillators
{
    enum class SHAPE {
        SINE,
        TRI,
        SAW,
        SQUARE,
        WHITE_NOISE,
        PINK_NOISE,
    };

    class SDSPOscillator
    {
    public:
        explicit SDSPOscillator(bool blep) : m_blep(blep) {

        }

        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate)
        {
            m_sampleRate = sampleRate;
            m_pinkingFilter.setCoefficients(SDSP::Filters::pinking().data());
        }

        float processSample() noexcept {
            float x = 0.0f;

            switch (m_currentShape) {
                case SHAPE::SINE: {
                   x = std::sin(m_phase * juce::MathConstants<float>::twoPi);
                   // No need to blep, sines don't alias..
                   break;
                }
                case SHAPE::TRI: {
                    x = std::asin(std::sin(m_phase* juce::MathConstants<float>::twoPi));
                    x /= juce::MathConstants<float>::halfPi;
                    // blep wise, this behaves as a square, but with a bit extra.. (which I dont fully get yet)
                    break;
                }
                case SHAPE::SAW: {
                    x = (2 * m_phase) - 1;
                    if(m_blep) {
                        x -= static_cast<float>(KMath::polyBlep(m_phase, m_phaseIncrement));
                    }
                    break;
                }
                case SHAPE::SQUARE: {
                    if(m_phase < 0.5f) x = -1.0f;
                    else if(m_phase > 0.5f) x = 1.0f;
                    if(m_blep) {
                        x += static_cast<float>(KMath::polyBlep(m_phase, m_phaseIncrement));
                        x -= static_cast<float>(KMath::polyBlep(std::fmod(m_phase + 0.5, 1), m_phaseIncrement));
                    }
                    break;
                }
                case SHAPE::WHITE_NOISE: {
                    x = juce::jmap<float>(juce::Random::getSystemRandom().nextFloat(), 0.0f, 1.0f, -1.0f, 1.0f);
                    break;
                }
                case SHAPE::PINK_NOISE: {
                    x = m_pinkingFilter.processSample(juce::jmap<float>(juce::Random::getSystemRandom().nextFloat(), 0.0f, 1.0f, -1.0f, 1.0f));
                    break;
                }
            }
            incrementPhase();
            return x;
        }

        [[maybe_unused]] SDSP_INLINE void setShape(SHAPE s) {
            m_currentShape = s;
            m_prevYn = 0.0f;
        }

        [[maybe_unused]] SDSP_INLINE void setFrequency(float newFrequency) noexcept {
            m_frequency = newFrequency;
            m_phaseIncrement = m_frequency / static_cast<float>(m_sampleRate);
        }

        [[maybe_unused]] SDSP_INLINE void retrigger() noexcept {
            m_phase = 0.0f;
        }

    private:
        SDSP_INLINE void incrementPhase() noexcept {
            m_phase += m_phaseIncrement;
            m_phase = std::fmod(m_phase, 1.0f);
        }

        bool m_blep{ true };
        double m_sampleRate{ 44100 };
        SHAPE m_currentShape{ SHAPE::SINE };
        float m_phase{ 0.0f };
        float m_phaseIncrement{ 0.0f };
        float m_frequency{ 0.0f };
        float m_prevYn{ 0.0f };

        BiquadCascade<1> m_pinkingFilter;
    };
}