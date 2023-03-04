#pragma once 
#include <juce_core/juce_core.h>
#include "../Macros.h"
namespace SDSP::Utils 
{ 
	//https://github.com/stekyne/PhaseVocoder/blob/master/DSP/Resample.h thank you again stekyne.
    SDSP_UNUSED inline static void linearResample(const float* const input, const int originalSize, float* const newSignal, const int newSignalSize)
    { 
	    const auto lerp = [&](float v0, float v1, float t)
	    {
	    	return (1.f - t) * v0 + t * v1;
	    };

	    // If the original signal is bigger than the new size, condense the signal to fit the new buffer
	    // otherwise expand the signal to fit the new buffer
	    const auto scale = static_cast<float>(originalSize) / static_cast<float>(newSignalSize);
	    float index = 0.f;

	    for (int i = 0; i < newSignalSize; ++i)
	    {
	    	const auto wholeIndex = (int)std::floor(index);
	    	const auto fractionIndex = index - static_cast<float>(wholeIndex);
	    	const auto sampleA = input[wholeIndex];
	    	const auto sampleB = input[wholeIndex + 1];
	    	newSignal[i] = lerp(sampleA, sampleB, fractionIndex);
	    	index += scale;
	    }
    }

    inline static juce::AudioBuffer<float> resampleAudioBuffer(juce::AudioBuffer<float>& toResample, double originalSampleRate, double newSampleRate) 
    { 
        juce::AudioBuffer<float> resampled;
        double ratio = originalSampleRate / newSampleRate;
        resampled.setSize(toResample.getNumChannels(), static_cast<int>(toResample.getNumSamples() / ratio));
        auto* read = toResample.getArrayOfReadPointers();
        auto* write = resampled.getArrayOfWritePointers();
        for (auto channel = 0; channel < toResample.getNumChannels(); channel++) {
            juce::LagrangeInterpolator resampler;
            resampler.reset();
            resampler.process(ratio, read[channel], write[channel], resampled.getNumSamples());
        }
        return resampled;
    }
}