#include <juce_dsp/juce_dsp.h> 
#include <functional>
namespace SDSP 
{ 
    class LFO 
    { 
    public: 
        enum class LFO_SHAPE { 
            SINE
        }
        
        LFO(LFO_SHAPE shape) { 
            m_oscillator.initialise(m_lookup[shape]);
        }
        
        ~LFO() { 

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
    private:
        bool m_hasBeenPrepared(false);
        float m_frequency{5.0f};
        juce::dsp::Oscillator<float> m_oscillator;
        static inline std::unordered_map<LFO_SHAPE, std::function<float(float)> > m_lookup { 
            {LFO_SHAPE::SINE, [](float x) { return juce::FastMathApproximations::sin<float>(x); }}
        };
    };
}