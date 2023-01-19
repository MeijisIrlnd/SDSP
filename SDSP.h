#pragma once
#if 0
BEGIN_JUCE_MODULE_DECLARATION

ID: SDSP
vendor: Syl
version: 1.0.0
name: SDSP
description: Syl DSP, catchy right?
license: Commercial
dependencies: juce_core, juce_dsp
END_JUCE_MODULE_DECLARATION
#endif
using APVTS = juce::AudioProcessorValueTreeState;
#include "Macros.h"
#include "Filters/Filters.h"
#include "Fourier/Fourier.h"
#include "Types/Types.h"
#include "Utils/Utils.h"
#include "Utils/LFO.h"
#include "KMath.h"
#include "Types/CircularBuffer.h"