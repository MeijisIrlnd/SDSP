#include "../Macros.h"
#include "../Oscillators/SDSPOscillator.h"

#include <juce_dsp/juce_dsp.h>

#include <functional>
#include <optional>
#include <unordered_map>

namespace SDSP 
{ 
    class LFO 
    { 
    public: 
        enum class LFO_SHAPE {
            SINE,
            PULSE,
        };
        
        explicit LFO(LFO_SHAPE shape) : m_oscillator(false) {
            m_oscillator.setShape(m_lookup[shape]);
        }

        explicit LFO(std::function<float(float)> shape) : m_oscillator(std::move(shape)) {

        }

        void prepareToPlay(int samplesPerBlockExpected, double sampleRate) { 
            m_oscillator.prepareToPlay(samplesPerBlockExpected, sampleRate);
            m_oscillator.setFrequency(m_frequency);
            m_hasBeenPrepared = true;
        }

        float processSample() {
            if (m_shapeToChange) {
                m_oscillator.setShape(m_lookup[m_shapeToChange.value()]);
                m_shapeToChange = std::nullopt;
            }

            if (m_functionToChange) {
                m_oscillator.setShape(std::move(m_functionToChange.value()));
                m_functionToChange = std::nullopt;
            }

            return m_oscillator.processSample();
        }

        void releaseResources() { 

        }

        void setFrequency(float newFrequency) { 
            m_frequency = newFrequency;
            if(m_hasBeenPrepared){ 
                m_oscillator.setFrequency(newFrequency);
            }
        }

        [[maybe_unused]] SDSP_INLINE void setShape(LFO_SHAPE newShape) {
            m_shapeToChange = newShape;
        }

        [[maybe_unused]] SDSP_INLINE void setShape(std::function<float(float)> newShape) {
            m_functionToChange = std::move(newShape);
        }

        SDSP_INLINE void retrigger() noexcept {
            m_oscillator.retrigger();
        }

        SDSP_INLINE void setPhase(float phase) noexcept {
            m_oscillator.setPhase(phase);
        }

    private:
        bool m_hasBeenPrepared{ false };
        float m_frequency{5.0f};
        Oscillators::SDSPOscillator m_oscillator;

        std::optional<LFO_SHAPE> m_shapeToChange = std::nullopt;
        std::optional<std::function<float(float)>> m_functionToChange = std::nullopt;

        std::unordered_map<LFO_SHAPE, std::function<float(float)>> m_lookup {
            {LFO_SHAPE::SINE, [](float x) { return juce::dsp::FastMathApproximations::sin<float>(x); }}
        };
    };
}