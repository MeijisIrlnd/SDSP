/*
  ==============================================================================

    CircularBuffer.h
    Created: 8 Jun 2021 1:24:59am
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <juce_core/juce_core.h>
#include "../Utils/LineRamp.h"
#include "../Macros.h"
namespace SDSP
{
    /// <summary>
    /// Mono circular buffer with a max size of 5 seconds, with tape delay style interpolation on delay time change
    /// </summary>
    /// <typeparam name="T"></typeparam>
    template <typename T>
    class CircularBuffer
    {
    public:

        CircularBuffer<T>(double interp, double sampleRate)
            : interpolationTimeMs(static_cast<int>(interp)),
            buffer(new T[static_cast<size_t>(sampleRate * 5)]),
            recordHead(buffer),
            playbackHead(buffer)
        {            
            std::memset(buffer, 0, static_cast<size_t>(sampleRate * 5) * sizeof(T));
        }

        explicit CircularBuffer<T>(double interp) : interpolationTimeMs(static_cast<int>(interp)) {

        }

        CircularBuffer<T>() : interpolationTimeMs(500) { }

        ~CircularBuffer()
        {
            if (buffer != nullptr) {
                delete[] buffer;
            }
        }
        
        // REFACTOR: what the fuck hahah
        void prepare(SDSP_UNUSED int spb, double sampleRate)
        {
            m_sampleRate = sampleRate;
            bufferSize = static_cast<size_t>(m_maxDelaySeconds * sampleRate);

            if (buffer != nullptr) { delete[] buffer; }
            buffer = new T[bufferSize];
            std::memset(buffer, 0, bufferSize * sizeof(T));
            
            recordHead = buffer;
            playbackHead = buffer;
            
            interpolator.prepare(sampleRate);
            interpolator.set(0, 0.5 * sampleRate, interpolationTimeMs);
        }

        void setMaxDelaySeconds(const double maxDelaySeconds) { 
            std::scoped_lock<std::mutex> sl(m_mutex);
            m_maxDelaySeconds = maxDelaySeconds;
            bufferSize = m_maxDelaySeconds * m_sampleRate;
            if (buffer != nullptr) { delete[] buffer; }
            buffer = new T[static_cast<size_t>(bufferSize)];
            std::memset(buffer, 0, bufferSize * sizeof(T));
            recordHead = buffer;
            playbackHead = buffer;
        }

        void clear()
        {
            std::memset(buffer, 0, bufferSize * sizeof(T));
            playbackHead = buffer;
            recordHead = buffer;

        }
        void setDelay(int m)
        {
            interpolator.set(m, interpolationTimeMs);
        }

        SDSP_UNUSED void setInterpolationRate(double interpTimeMs)
        {
            interpolationTimeMs = static_cast<int>(interpTimeMs);
        }

        inline T getNextSample(T in)
        {
            std::scoped_lock<std::mutex> sl(m_mutex);
            recordHead = playbackHead - (int)interpolator.process();
            while (recordHead < buffer) {
                recordHead += bufferSize;
            }
            
            T y;
            *playbackHead = in;
            ++playbackHead;
            y = *recordHead;
            ++recordHead;
            if ((playbackHead - buffer) >= bufferSize) playbackHead -= bufferSize;
            if ((recordHead - buffer) >= bufferSize) recordHead -= bufferSize;
            return y;

        }

    private:
        std::mutex m_mutex;
        juce::SmoothedValue<float> delayTimeInterpolator;
        LineRamp<double> interpolator;
        T* buffer = nullptr;
        T* recordHead;
        T* playbackHead;
        double m_sampleRate{};
        size_t bufferSize{};
        int interpolationTimeMs = 500;
        double m_maxDelaySeconds{5};

    };
}