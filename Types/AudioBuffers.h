//
// Created by Syl on 12/04/2023.
//
#pragma once
#include <juce_audio_basics/juce_audio_basics.h>
namespace SDSP
{
    struct BufferWithSampleRate
    {
        BufferWithSampleRate() = default;
        BufferWithSampleRate(juce::AudioBuffer<float>&& buffer, double sampleRate) :
                m_buffer(buffer), m_sampleRate(sampleRate) {}
        juce::AudioBuffer<float> m_buffer;
        double m_sampleRate = 0.0;
    };

    class BufferTransfer
    {
    public:
        void set(BufferWithSampleRate&& buffer) {
            const juce::SpinLock::ScopedLockType lock(m_mutex);
            m_buffer = std::move(buffer);
            m_newBuffer = true;
        }

        template<typename Fn>
        void get(Fn&& fn) {
            const juce::SpinLock::ScopedTryLockType lock(m_mutex);
            if (lock.isLocked() && m_newBuffer) {
                fn(m_buffer);
                m_newBuffer = false;
            }
        }

    private:
        BufferWithSampleRate m_buffer;
        bool m_newBuffer = false;
        juce::SpinLock m_mutex;
    };

}