#pragma once 
#include <juce_dsp/juce_dsp.h>
#include "../Types/BlockCircularBuffer.h"
#include "../Utils/Resampling.h"
namespace SDSP::Fourier 
{ 
    /// <summary>
    /// More lifting from stekyne, https://github.com/stekyne/PhaseVocoder/blob/master/DSP/PhaseVocoder.h for the stft process snippet..
    /// Window length should be a POT, 
    /// </summary>
    
    template<int WindowSize, int FFTOrder>
    class STFT
    {
    public:
        STFT() : m_fft(FFTOrder), 
            m_analysisBuffer(WindowSize), 
            m_synthesisBuffer(WindowSize * 3), 
            m_windowData(WindowSize)
        {
            juce::dsp::WindowingFunction<float>::fillWindowingTables(m_windowData.data(), WindowSize, juce::dsp::WindowingFunction<float>::hann, false);
            m_spectralBuffer.resize(m_spectralBufferSize);
            std::fill(m_spectralBuffer.begin(), m_spectralBuffer.end(), 0.0f);
            
        }
        SDSP_INLINE constexpr static int getFFTSize() { return m_fftSize; }
        SDSP_INLINE constexpr static int getAnalysisHopSize() { return m_analysisHopSize; }
        SDSP_INLINE constexpr static int getWindowSize() { return WindowSize; }

        SDSP_INLINE void process(const float* data, float* out, size_t bufferSize, std::function<void(float*, size_t)>& callback) noexcept
        {
            for (auto internalOffset = 0, internalBufferSize = 0; internalOffset < bufferSize; internalOffset += internalBufferSize)
            {
                const auto remainingSamples = (bufferSize - internalOffset);
                internalBufferSize = m_incomingSampleCount + remainingSamples >= m_samplesUntilNextProcess ? m_samplesUntilNextProcess - m_incomingSampleCount : remainingSamples;
                const auto previousAnalysisWriteIndex = m_analysisBuffer.getReadIndex();
                m_analysisBuffer.write(data + internalOffset, internalBufferSize);
                m_incomingSampleCount += internalBufferSize;
                if (m_incomingSampleCount >= m_samplesUntilNextProcess)
                {
                    m_isProcessing = true;
                    m_incomingSampleCount -= m_samplesUntilNextProcess;
                    m_samplesUntilNextProcess = m_analysisHopSize;
                    auto spectralBufferData = m_spectralBuffer.data();
                    m_analysisBuffer.setReadHopSize(m_analysisHopSize);
                    m_analysisBuffer.read(spectralBufferData, WindowSize);
                    juce::FloatVectorOperations::multiply(spectralBufferData, m_windowData.data(), WindowSize);
                    std::rotate(spectralBufferData, spectralBufferData + (WindowSize / 2), spectralBufferData + WindowSize);
                    m_fft.performRealOnlyForwardTransform(spectralBufferData);
                    callback(spectralBufferData, m_spectralBufferSize);
                    m_fft.performRealOnlyInverseTransform(spectralBufferData);
                    std::rotate(spectralBufferData, spectralBufferData + (WindowSize / 2), spectralBufferData + WindowSize);
                    juce::FloatVectorOperations::multiply(spectralBufferData, m_windowData.data(), WindowSize);
                    m_synthesisBuffer.setWriteHopSize(m_synthesisHopSize);
                    m_synthesisBuffer.overlapWrite(spectralBufferData, WindowSize);
                }
                if (!m_isProcessing)
                {
                    std::fill(out + internalOffset, out + internalOffset + internalBufferSize, 0.0f);
                    continue;
                }
                m_synthesisBuffer.read(out + internalOffset, internalBufferSize);
            }
        }



    private:
        juce::dsp::FFT m_fft;
        constexpr static int m_fftSize = 1 << FFTOrder;

        BlockCircularBuffer<float> m_analysisBuffer, m_synthesisBuffer;
        int m_incomingSampleCount{ 0 };
        int m_samplesUntilNextProcess{WindowSize};
        bool m_isProcessing{ false };

        constexpr static int m_analysisHopSize = m_fftSize / 4;
        constexpr static int m_synthesisHopSize = m_fftSize / 4;
        constexpr static int m_spectralBufferSize = WindowSize * 2;

        std::vector<float> m_spectralBuffer;
        std::vector<float> m_windowData;
    };

    template<int WindowSize, int FFTOrder> 
    class STFTWrapper 
    { 
    public:
        STFTWrapper() 
        {
        }
        ~STFTWrapper() 
        { 
        }

        STFTWrapper(const STFTWrapper& other) = delete;
        STFTWrapper(STFTWrapper&& other) noexcept = delete;
        STFTWrapper& operator=(const STFTWrapper& other) = delete;
        STFTWrapper& operator=(STFTWrapper&& other) noexcept = delete;

        float processSample(float in, const size_t bufferSize) 
        { 
            if(m_data.size() != bufferSize) { 
                m_data.resize(bufferSize);
                juce::FloatVectorOperations::fill(m_data.data(), 0.0f, bufferSize);
            }
            if(m_accumulator.size() != bufferSize) { 
                m_accumulator.resize(bufferSize);
                juce::FloatVectorOperations::fill(m_accumulator.data(), 0.0f, bufferSize);
            }
            if(!m_hasPerformedFirstTransform) { 
                m_samplesUntilUpdate = bufferSize;
                m_hasPerformedFirstTransform = true;
            }
            if(m_samplesUntilUpdate == 0) { 
                m_stft.process(m_accumulator.data(), m_data.data(), bufferSize, stftCallback);
                m_samplesUntilUpdate = bufferSize;
                m_writePos = 0;
            }
            m_accumulator[m_writePos] = in;
            float out = m_data[m_writePos];
            ++m_writePos;
            --m_samplesUntilUpdate;
            return out;
        }

        std::function<void(float*, size_t)> stftCallback{nullptr};
    private: 
        bool m_hasPerformedFirstTransform{false};
        std::vector<float> m_accumulator, m_data;

        int m_samplesUntilUpdate{0}, m_writePos{0};
        STFT<WindowSize, FFTOrder> m_stft;
    };
}