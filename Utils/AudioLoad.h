#pragma once 
#include <juce_core/juce_core.h>
#include "Resampling.h"
namespace SDSP::Utils 
{ 
    static inline double loadAudioFile(juce::AudioFormatManager& formatManager, const juce::File& toLoad, juce::AudioBuffer<float>& destBuffer) 
    { 
        std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(toLoad));
        destBuffer.setSize(reader->lengthInSamples, static_cast<int>(reader->lengthInSamples));
        reader->read(&destBuffer, 0, static_cast<int>(reader->lengthInSamples), 0, true, true);
        return reader->sampleRate;
    }

    static inline void loadAudioFileWithResample(juce::AudioFormatManager& formatManager, const juce::File& toLoad, const double sampleRate, juce::AudioBuffer<float>& destBuffer)
    {
        juce::AudioBuffer<float> loadBuffer;
        auto loadedSampleRate = loadAudioFile(formatManager, toLoad, loadBuffer);
        destBuffer = resampleAudioBuffer(loadBuffer, loadedSampleRate, sampleRate);
    }

    static inline double loadFromMemory(juce::AudioFormatManager& formatManager, const void* data, int size, juce::AudioBuffer<float>& destBuffer) 
    { 
        std::unique_ptr<juce::MemoryInputStream> ip(new juce::MemoryInputStream(data, size, false));
        std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(std::move(ip)));
        if (reader.get() == nullptr) return -1;
        destBuffer.setSize(reader->numChannels, static_cast<int>(reader->lengthInSamples), false);
        reader->read(&destBuffer, 0, static_cast<int>(reader->lengthInSamples), 0, true, true);
        return reader->sampleRate;
    }

    static inline void loadFromMemoryWithResample(juce::AudioFormatManager& formatManager, const void* data, int size, const double sampleRate, juce::AudioBuffer<float>& destBuffer)
    {
        juce::AudioBuffer<float> loadBuffer;
        auto loadedSampleRate = loadFromMemory(formatManager, data, size, loadBuffer);
        destBuffer = resampleAudioBuffer(loadBuffer, loadedSampleRate, sampleRate);
    }
}