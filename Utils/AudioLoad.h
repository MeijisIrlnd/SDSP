#pragma once 
#include <juce_core/juce_core.h>
#include "Resampling.h"
namespace SDSP::Utils 
{ 
    static inline double loadAudioFile(juce::AudioFormatManager& formatManager, const juce::File& toLoad, juce::AudioBuffer<float>& destBuffer) 
    { 
        std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(toLoad));
        destBuffer.setSize(2, static_cast<int>(reader->lengthInSamples));
        reader->read(&destBuffer, 0, static_cast<int>(reader->lengthInSamples), 0, true, true);
        return reader->sampleRate;
    }

    static inline void loadAudioFileWithResample(juce::AudioFormatManager& formatManager, const juce::File& toLoad, const double sampleRate, juce::AudioBuffer<float>& destBuffer)
    {
        juce::AudioBuffer<float> loadBuffer;
        auto loadedSampleRate = loadAudioFile(formatManager, toLoad, loadBuffer);
        destBuffer = resampleAudioBuffer(loadBuffer, loadedSampleRate, sampleRate);
    }
}