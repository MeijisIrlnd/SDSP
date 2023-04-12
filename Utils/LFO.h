#include <juce_dsp/juce_dsp.h> 
#include <unordered_map>
#include <functional>
namespace SDSP 
{ 
    class LFO 
    { 
    public: 
        enum class LFO_SHAPE {
            SINE
        };
        
        explicit LFO(const LFO_SHAPE& shape) {
            m_oscillator.initialise(m_lookup[shape]);
        }

        explicit LFO(const std::function<float(float)>& shape) {
            m_oscillator.initialise(shape);
        }



        void prepareToPlay(int samplesPerBlockExpected, double sampleRate) { 
            m_oscillator.prepare({sampleRate, static_cast<juce::uint32>(samplesPerBlockExpected), 2});
            m_oscillator.setFrequency(m_frequency);
            m_hasBeenPrepared = true;
        }

        float processSample() { 
            return m_oscillator.processSample(0.0f);
        }

        void releaseResources() { 
            m_oscillator.reset();
        }

        void setFrequency(float newFrequency) { 
            m_frequency = newFrequency;
            if(m_hasBeenPrepared){ 
                m_oscillator.setFrequency(newFrequency);
            }
        }

        [[maybe_unused]] SDSP_INLINE void setShape(LFO_SHAPE newShape) {
            m_oscillator.initialise(m_lookup[newShape]);
        }

        [[maybe_unused]] SDSP_INLINE void setShape(const std::function<float(float)>& newShape) {
            m_oscillator.initialise(newShape);
        }

    private:
        bool m_hasBeenPrepared{ false };
        float m_frequency{5.0f};
        juce::dsp::Oscillator<float> m_oscillator;
        std::unordered_map<LFO_SHAPE, std::function<float(float)> > m_lookup { 
            {LFO_SHAPE::SINE, [](float x) { return juce::dsp::FastMathApproximations::sin<float>(x); }}
        };
    };
}