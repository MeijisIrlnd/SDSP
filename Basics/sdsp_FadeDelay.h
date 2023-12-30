#ifndef SDSP_FADEDELAY_H
#define SDSP_FADEDELAY_H
#include <juce_dsp/juce_dsp.h>
namespace SDSP::Basics {
    class FadeDelay {
    public:
        explicit FadeDelay(int xfadeTimeMs) : m_xfadeSeconds(static_cast<float>(xfadeTimeMs) / 1000.0f) {
            //
        }

        void prepare(const juce::dsp::ProcessSpec& spec) {
            for (auto& d : m_delayLines) {
                d.prepare(spec);
            }
            const auto xFadeTimeSamples = static_cast<int>(m_xfadeSeconds * static_cast<float>(spec.sampleRate));
            m_indexInterpolation.reset(xFadeTimeSamples);
            m_indexInterpolation.setCurrentAndTargetValue(0.0f);
        }

        void pushSample(int channel, float x) noexcept {
            m_delayLines[0].pushSample(channel, x);
            m_delayLines[1].pushSample(channel, x);
        }

        [[nodiscard]] float popSample(int channel) noexcept {
            const auto xFade = m_indexInterpolation.getNextValue();
            const auto outA = m_delayLines[0].popSample(channel);
            const auto outB = m_delayLines[1].popSample(channel);
            return (outA * (1.0f - xFade)) + (outB * xFade);
        }

        [[nodiscard]] float processSample(int channel, float x) noexcept {
            const auto out = popSample(channel);
            pushSample(channel, x);
            return out;
        }

        void reset() noexcept {
            for (auto& d : m_delayLines) {
                d.reset();
            }
        }

        void setMaximumDelayInSamples(int newMaxSamples) noexcept {
            for (auto& d : m_delayLines) {
                d.setMaximumDelayInSamples(newMaxSamples);
            }
        }

        void setDelay(float newDelaySamples) noexcept {
            if (m_delayLines[m_index].getDelay() == newDelaySamples) return;
            if (m_indexInterpolation.isSmoothing()) return;
            if (m_firstSet) {
                m_delayLines[m_index].setDelay(newDelaySamples);
                m_firstSet = false;
                return;
            }
            ++m_index;
            m_index %= 2;
            m_delayLines[m_index].setDelay(newDelaySamples);
            m_indexInterpolation.setTargetValue(static_cast<float>(m_index));
        }

    private:
        std::array<juce::dsp::DelayLine<float, juce::dsp::DelayLineInterpolationTypes::None>, 2> m_delayLines;
        // blend between these delay lines on a set..
        size_t m_index{ 0 };
        juce::SmoothedValue<float> m_indexInterpolation{};
        bool m_firstSet{ true };
        const float m_xfadeSeconds;
    };
} // namespace SDSP::Basics
#endif