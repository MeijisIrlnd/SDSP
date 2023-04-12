//
// Created by Syl on 12/04/2023.
//
#pragma once
#include "CircularBuffer.h"
#include <array>
namespace SDSP
{
    class CascadingCircularBuffer
    {
    public:
        CascadingCircularBuffer()
        {

        }

        CascadingCircularBuffer(const CascadingCircularBuffer& other) = default;

        void prepare(int samplesPerBlockExpected, double sampleRate)
        {
            for (auto i = 0; i < 8; i++) {
                m_buffers[i].prepare(samplesPerBlockExpected, sampleRate);

            }
        }

        float getNextSample(float in)
        {
            float y = 0;
            for (auto i = 0; i < 8; i++) {
                y += m_buffers[i].getNextSample(in);
            }
            return (y / 8.0f);
        }

        void getNextSample(float in, std::array<float, 8>& dest)
        {
            for (auto tap = 0; tap < 8; tap++) {
                dest[tap] = m_buffers[tap].getNextSample(in);
            }
        }

        void setDelayTime(int timeSamples)
        {
            for (auto i = 0; i < 8; i++) {
                m_buffers[i].setDelay(static_cast<int>(timeSamples / std::pow(2, i)));
            }
        }
    private:

        float test = 0;
        std::array<CircularBuffer<float>, 8> m_buffers;
    };
}