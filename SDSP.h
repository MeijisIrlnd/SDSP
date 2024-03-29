#pragma once
#if 0
BEGIN_JUCE_MODULE_DECLARATION

ID: SDSP
vendor: Syl
version: 1.0.0
name: SDSP
description: Syl DSP, catchy right?
license: Commercial
dependencies: juce_audio_utils, juce_core, juce_dsp
END_JUCE_MODULE_DECLARATION
#endif

#include "Macros.h"
#include "Filters/Filters.h"
#include "Utils/Utils.h"
#include "KMath.h"
#include "Types/CircularBuffer.h"
#include "Types/MultitapCircularBuffer.h"
#include "Helpers/Helpers.h"
#include "Utils/DelayTimeInterpolator.h"
#include "Utils/BlendableLFO.h"
#include "Oscillators/SDSPOscillator.h"