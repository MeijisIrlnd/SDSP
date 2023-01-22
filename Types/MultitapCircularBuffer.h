#include <juce_audio_basics/juce_audio_basics.h>
// TODO: TEST AND MODULATE
namespace SDSP 
{ 
    template<int NTAPS>
    class MultitapCircularBuffer 
    { 
        public: 
            MultitapCircularBuffer() { 
                std::fill(m_delayTimes.begin(), m_delayTimes.end(), 0.0f);
            }
            void prepareToPlay(int samplesPerBlockExpected, double sampleRate) { 
                m_sampleRate = sampleRate;
                m_buffer.resize(m_maxDelaySeconds * sampleRate);
                std::fill(m_buffer.begin(), m_buffer.end(), 0.0f);
                std::fill(m_read.begin(), m_read.end(), 0.0f);
                m_hasBeenPrepared = true;
            }

            std::array<float, NTAPS> processSample(float x) { 
                juce::FloatVectorOperations::fill(m_output.data(), 0.0f, NTAPS);
                m_buffer[m_write] = x;
                for(auto tap = 0; tap < NTAPS; tap++) { 
                    m_read[tap] = m_write - m_delayTimes[tap];
                    while(m_read[tap] < m_buffer.size()) { 
                        m_read[tap] += m_buffer.size();
                    }
                    m_output[tap] = m_buffer[m_read[tap]];
                    ++m_read[tap];
                    if(m_read[tap] >= m_buffer.size() - 1) { 
                        m_read[tap] -= m_buffer.size();
                    }
                }
                ++m_write;
                if(m_write >= m_buffer.size() - 1) { 
                    m_write -= m_buffer.size();
                }
                return m_output;
            }

            float releaseResources { 

            }

            void setMaximumDelaySeconds(float maxDelaySeconds) { 
                m_maxDelaySeconds = m_maxDelaySeconds;
                if(m_hasBeenPrepared) { 
                    m_buffer.resize(m_maxDelaySeconds * m_sampleRate);
                    std::fill(m_buffer.begin(), m_buffer.end(), 0.0f);
                }
            }

            void setDelayTime(int delayTime, int tap) { 
                m_delayTimes[tap] = delayTime;
            }

        private:
            bool m_hasBeenPrepared{false};
            double m_sampleRate;
            std::vector<float> m_buffer;
            size_t m_write = 0;
            std::array<size_t, NTAPS> m_read;
            std::array<int, NTAPS> m_delayTimes;
            float m_maxDelaySeconds{2};
            std::array<float, NTAPS> m_output;
    };      
}