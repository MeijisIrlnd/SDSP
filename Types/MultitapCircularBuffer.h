#include <juce_audio_basics/juce_audio_basics.h>
// TODO: TEST AND MODULATE
namespace SDSP {
    template <int NTAPS>
    class [[maybe_unused]] MultitapCircularBuffer {
    public:
        MultitapCircularBuffer() {
            std::fill(m_delayTimes.begin(), m_delayTimes.end(), 0);
        }
        void prepareToPlay(SDSP_UNUSED int samplesPerBlockExpected, double sampleRate) {
            m_sampleRate = sampleRate;
            m_buffer.resize(static_cast<size_t>(static_cast<double>(m_maxDelaySeconds) * sampleRate));
            std::fill(m_buffer.begin(), m_buffer.end(), 0.0f);
            std::fill(m_read.begin(), m_read.end(), 0);
            for (size_t tap = 0; tap < static_cast<size_t>(NTAPS); tap++) {
                m_smoothedDelayTimes[tap].reset(sampleRate, 0.5f);
                m_smoothedDelayTimes[tap].setCurrentAndTargetValue(static_cast<float>(m_delayTimes[tap]));
            }
            m_hasBeenPrepared = true;
        }

        std::array<float, NTAPS> processSample(float x) {
            juce::FloatVectorOperations::fill(m_output.data(), 0.0f, NTAPS);
            while (m_write >= m_buffer.size()) {
                m_write -= static_cast<int>(m_buffer.size());
            }

            m_buffer[m_write] = x;
            for (size_t tap = 0; tap < static_cast<size_t>(NTAPS); ++tap) {
                m_read[tap] = m_write - static_cast<int>(m_smoothedDelayTimes[tap].getNextValue());
                // m_read[tap] = m_write - m_delayTimes[tap];
                while (m_read[tap] < 0) {
                    m_read[tap] += static_cast<int>(m_buffer.size());
                }
                // while(m_read[tap] < m_buffer.size()) {
                //     m_read[tap] += m_buffer.size();
                // }
                m_output[tap] = m_buffer[m_read[tap]];
                ++m_read[tap];
                if (m_read[tap] >= static_cast<int>(m_buffer.size()) - 1) {
                    m_read[tap] -= static_cast<int>(m_buffer.size());
                }
            }
            ++m_write;

            return m_output;
        }

        float releaseResources() {
        }

        SDSP_UNUSED void setMaximumDelaySeconds(float maxDelaySeconds) {
            m_maxDelaySeconds = maxDelaySeconds;
            if (m_hasBeenPrepared) {
                m_buffer.resize(m_maxDelaySeconds * m_sampleRate);
                std::fill(m_buffer.begin(), m_buffer.end(), 0.0f);
            }
        }

        void setDelayTime(int delayTime, int tap) {
            m_delayTimes[tap] = delayTime;
            if (m_hasBeenPrepared) {
                m_smoothedDelayTimes[tap].setTargetValue(static_cast<float>(delayTime));
            }
        }

        void clearBuffers() {
            std::fill(m_buffer.begin(), m_buffer.end(), 0.0f);
        }

    private:
        bool m_hasBeenPrepared{ false };
        double m_sampleRate{};
        std::vector<float> m_buffer;
        int m_write = 0;
        std::array<int, NTAPS> m_read;
        std::array<int, NTAPS> m_delayTimes;
        float m_maxDelaySeconds{ 2 };
        std::array<float, NTAPS> m_output;
        std::array<juce::SmoothedValue<float>, NTAPS> m_smoothedDelayTimes;
    };
} // namespace SDSP