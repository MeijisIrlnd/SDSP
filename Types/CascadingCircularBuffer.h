//
// Created by Syl on 12/04/2023.
//
#pragma once
#include "CircularBuffer.h"
#include <array>
namespace SDSP
{
    class [[maybe_unused]] CascadingCircularBuffer
    {
    public:
        CascadingCircularBuffer() = default;

        CascadingCircularBuffer(const CascadingCircularBuffer& other) = default;

        void prepare(int samplesPerBlockExpected, double sampleRate)
        {
            for (size_t i = 0; i < 8; i++) {
                m_buffers[i].prepare(samplesPerBlockExpected, sampleRate);

            }
        }

        float getNextSample(float in)
        {
            float y = 0;
            for (size_t i = 0; i < 8; i++) {
                y += m_buffers[i].getNextSample(in);
            }
            return (y / 8.0f);
        }

        void getNextSample(float in, std::array<float, 8>& dest)
        {
            for (size_t tap = 0; tap < 8; tap++) {
                dest[tap] = m_buffers[tap].getNextSample(in);
            }
        }

        void setDelayTime(int timeSamples)
        {
            for (auto i = 0; i < 8; i++) {
                m_buffers[static_cast<size_t>(i)].setDelay(static_cast<int>(timeSamples / std::pow(2, i)));
            }
        }

        void clearBuffers() {
            for(auto& b : m_buffers) {
                b.clear();
            }
        }
    private:
        std::array<CircularBuffer<float>, 8> m_buffers;
    };
}