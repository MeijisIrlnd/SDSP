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
        SAW,
        SQUARE,
        TRI,
        PULSE,
        WHITE_NOISE,
        PINK_NOISE,
        CUSTOM,
    };

    class [[maybe_unused]] SDSPOscillator
    {
    public:

        using OscFunction = std::function<float(float)>;

        [[maybe_unused]] explicit SDSPOscillator(bool blep) : m_blep(blep) {

        }

        [[maybe_unused]] explicit SDSPOscillator(OscFunction function) : m_blep(false), m_customFunction(std::move(function)) {

        }

        [[maybe_unused]] SDSPOscillator(bool blep, SHAPE shape) : m_blep(blep), m_currentShape(shape) {

        }

        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate) {
            m_sampleRate = sampleRate;
            m_pinkingFilter.setCoefficients(SDSP::Filters::pinking().data());
        }

        float processSample() noexcept {
            if(m_phaseIncrement == 0.0f) return 0.0f;
            float x = 0.0f;
            float offsetPhase = std::fmod(m_phase + m_offset, 1.0f);

            switch (m_currentShape) {
                case SHAPE::SINE: {
                   x = std::sin(offsetPhase * juce::MathConstants<float>::twoPi);
                   // No need to blep, sines don't alias..
                   break;
                }
                case SHAPE::TRI: {
                    x = std::asin(std::sin(offsetPhase * juce::MathConstants<float>::twoPi));
                    x /= juce::MathConstants<float>::halfPi;
                    // blep wise, this behaves as a square, but with a bit extra.. (which I dont fully get yet)
                    break;
                }
                case SHAPE::SAW: {
                    x = (2 * offsetPhase) - 1;
                    if(m_blep) {
                        x -= static_cast<float>(KMath::polyBlep(offsetPhase, m_phaseIncrement));
                    }
                    break;
                }
                case SHAPE::SQUARE: {
                    if(offsetPhase < 0.5f) x = -1.0f;
                    else if(offsetPhase > 0.5f) x = 1.0f;
                    if(m_blep) {
                        x += static_cast<float>(KMath::polyBlep(offsetPhase, m_phaseIncrement));
                        x -= static_cast<float>(KMath::polyBlep(std::fmod(offsetPhase + 0.5, 1), m_phaseIncrement));
                    }
                    break;
                }
                case SHAPE::PULSE: {
                    if (offsetPhase < m_pulseWidth) {
                        x = -1.0f;
                    } else {
                        x = 1.0f;
                    }

                    if (m_blep) {
                        x += static_cast<float>(KMath::polyBlep(offsetPhase, m_phaseIncrement));
                        x -= static_cast<float>(KMath::polyBlep(std::fmod(offsetPhase + 0.5, 1), m_phaseIncrement));
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
                case SHAPE::CUSTOM: {
                    if (m_tempFunction) {
                        m_customFunction = std::move(*m_tempFunction);
                        m_tempFunction.reset();
                    }
                    x = m_customFunction(offsetPhase);
                    break;
                }
            }
            incrementPhase();
            return x;
        }

        [[maybe_unused]] SDSP_INLINE void setShape(SHAPE s) {
            m_currentShape = s;
        }

        [[maybe_unused]] [[nodiscard]] SDSP_INLINE SHAPE getShape() const { return m_currentShape; }
        [[maybe_unused]] SDSP_INLINE void setFunction(OscFunction function) {
            m_tempFunction = std::move(function);
            m_currentShape = SHAPE::CUSTOM;
        }

        [[maybe_unused]] SDSP_INLINE void setFrequency(float newFrequency) noexcept {
            m_frequency = newFrequency;
            m_phaseIncrement = m_frequency / static_cast<float>(m_sampleRate);
        }

        [[maybe_unused]] SDSP_INLINE void retrigger() noexcept {
            m_phase = 0.0f;
        }

        [[maybe_unused]] SDSP_INLINE void setPulsewidth(float pulsewidth) noexcept {
            m_pulseWidth = pulsewidth;
        }

        [[maybe_unused]] SDSP_INLINE void setPhase(float phase) noexcept {
            m_phase = phase;
        }

        [[maybe_unused]] [[nodiscard]] SDSP_INLINE float getPhase() const noexcept {
            return m_phase;
        }

        // incoming value expect to be in range 0 - pi
        [[maybe_unused]] SDSP_INLINE void setOffset(float offset) noexcept {
            m_offset = juce::jmap(offset, 0.0f, juce::MathConstants<float>::pi, 0.0f, 1.0f);
        }

        // DANGER
        SDSP_INLINE void setPhaseIncrement(float newPhaseIncrement) {
            m_phaseIncrement = newPhaseIncrement;
        }

    private:

        SDSP_INLINE void incrementPhase() noexcept {
            m_phase += m_phaseIncrement;
            m_phase = std::fmod(m_phase, 1.0f);
        }

        SHAPE m_currentShape{ SHAPE::SINE };
        OscFunction m_customFunction;
        std::optional<OscFunction> m_tempFunction{ std::nullopt };

        bool m_blep{ true };
        double m_sampleRate{ 44100 };
        float m_phase{ 0.0f };
        float m_phaseIncrement{ 0.0f };
        float m_offset{ 0.0f };
        float m_frequency{ 0.0f };
        float m_pulseWidth{ 0.5f };

        BiquadCascade<1> m_pinkingFilter;
    };
}