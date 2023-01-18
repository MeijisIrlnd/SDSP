/*
  ==============================================================================

    LineRamp.h
    Created: 8 Jun 2021 1:25:58am
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <juce_core/juce_core.h>
namespace SDSP
{
    /// <summary>
    /// Linear interpolator from point start to point end
    /// </summary>
    /// <typeparam name="T"></typeparam>
    template <typename T>
    class LineRamp
    {
    public:
        LineRamp<T>()
        {

        }
        virtual ~LineRamp() {}
        void prepare(double sampleRate)
        {
            samplesPerMs = sampleRate / 1000.0;
        }

        virtual T process()
        {
            T output = currentValue;
            if (remainingTicks > 0)
            {
                currentValue += increment;
                --remainingTicks;
                if (remainingTicks <= 0)
                {
                    currentValue = targetValue;
                    running = false;
                }
            }
            return output;
        }

        void set(double timeMs)
        {
            set(startValue, targetValue, timeMs);
        }

        void set(T target, double timeMs)
        {
            startValue = currentValue;
            set(startValue, target, timeMs);
        }

        void set(T start, T target, double timeMs)
        {
            startValue = start;
            targetValue = target;
            totalTicks = static_cast<juce::uint64>(samplesPerMs * timeMs);
            remainingTicks = totalTicks;
            currentValue = startValue;
            increment = (targetValue - startValue) / (double)totalTicks;
            running = true;
        }

        void stop()
        {
            remainingTicks = 0;
            running = false;
        }

    protected:
        bool running = false;
        T currentValue = 0.0;
        double samplesPerMs;
        juce::uint64 remainingTicks;
        juce::uint64 totalTicks;
        double increment;
        T targetValue;
        T startValue;
    };
}
