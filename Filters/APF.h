//
// Created by Syl on 05/03/2023.
//

#pragma once
#include "../Oscillators/SDSPOscillator.h"
#include <juce_dsp/juce_dsp.h>
namespace SDSP::Filters
{
    class [[maybe_unused]] APFBase {
    public:
        [[maybe_unused]] float tap(const int index) {
            return m_delayLine.popSample(0, static_cast<float>(index), false);
        }
        [[maybe_unused]] virtual void setDelayTimeSamples(int newDelayTimeSamples) {
            m_delayTimeSamples = newDelayTimeSamples;
            if(m_hasBeenPrepared) {
                m_delayLine.setDelay(static_cast<float>(m_delayTimeSamples));
            }
        }

        [[maybe_unused]] void setCoeff(float newCoeff) {
            m_coeff = newCoeff;
            if(m_hasBeenPrepared) {
                m_smoothedCoeff.setTargetValue(m_coeff);
            }
        }

        virtual void prepareToPlay(int samplesPerBlockExpected, double sampleRate) {
            m_sampleRate = sampleRate;
            m_hasBeenPrepared = true;
            juce::dsp::ProcessSpec spec{sampleRate, static_cast<juce::uint32>(samplesPerBlockExpected), 1};
            m_delayLine.prepare(spec);
            // Allow pad for delay modulation
            m_delayLine.setMaximumDelayInSamples(static_cast<int>(sampleRate / 4.0));
            m_delayLine.setDelay(static_cast<float>(m_delayTimeSamples));
            m_smoothedCoeff.reset(sampleRate, 0.1);
            m_smoothedCoeff.setCurrentAndTargetValue(m_coeff);
        }

        virtual float processSample(float in) = 0;

        void releaseResources()
        {

        }

        virtual void clearBuffers() noexcept = 0;

    protected:
        double m_sampleRate{ 44100 };
        bool m_hasBeenPrepared{ false };
        int m_delayTimeSamples{ 0 };
        juce::dsp::DelayLine<float> m_delayLine;
        float m_coeff{ 0.5f };
        juce::SmoothedValue<float> m_smoothedCoeff;
    };

    class [[maybe_unused]] APF : public APFBase {
    public:
        float processSample(float in) override
        {
            auto delayed = m_delayLine.popSample(0);
            in += (delayed * -m_coeff);
            m_delayLine.pushSample(0, in);
            auto feedforward = in * m_coeff;
            return (delayed + feedforward);
        }

        void clearBuffers() noexcept override {
            m_delayLine.reset();
        }
    };

    class [[maybe_unused]] ModulatedAPF : public APFBase {
    public:
        ModulatedAPF() : m_lfo(false) {
            m_lfo.setShape(Oscillators::SHAPE::SINE);
        }

        [[maybe_unused]] void setLfoRate(const float newRate) noexcept {
            m_lfo.setFrequency(newRate);
        }

        [[maybe_unused]] void setExcursionSeconds(float newExcursionSeconds) {
            m_excursionSeconds = newExcursionSeconds;
        }

        void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override
        {
            APFBase::prepareToPlay(samplesPerBlockExpected, sampleRate);
            juce::dsp::ProcessSpec spec{sampleRate / static_cast<double>(m_updateRate), static_cast<juce::uint32>(samplesPerBlockExpected), 1};
            m_lfo.prepareToPlay(samplesPerBlockExpected, sampleRate);
            m_smoothedDelayTimeSamples.reset(sampleRate, 0.1);
            m_smoothedDelayTimeSamples.setCurrentAndTargetValue(static_cast<float>(m_delayTimeSamples));
        }

        float processSample(float in) override {
            if(m_samplesUntilUpdate == 0) {
                auto lfoSample = m_lfo.processSample();
                auto excursionSamples = static_cast<float>(m_excursionSeconds * m_sampleRate);
                lfoSample = juce::jmap<float>(lfoSample, -1, 1, -excursionSamples, excursionSamples);
                auto newDelayTime = (static_cast<float>(m_delayTimeSamples) + lfoSample);
                m_smoothedDelayTimeSamples.setTargetValue(newDelayTime);
                m_samplesUntilUpdate = m_updateRate;
            }
            --m_samplesUntilUpdate;
            // in case you need to revert, here's the inverted coeff version
//            auto currentDelay = m_smoothedDelayTimeSamples.getNextValue();
//            auto delayed = m_delayLine.popSample(0, currentDelay);
//            in += (delayed * m_coeff);
//            m_delayLine.pushSample(0, in);
//            auto feedforward = in * -m_coeff;
//            return (delayed + feedforward);

            auto currentDelay = m_smoothedDelayTimeSamples.getNextValue();
            auto delayed = m_delayLine.popSample(0, currentDelay);
            in += (delayed * -m_coeff);
            m_delayLine.pushSample(0, in);
            auto feedforward = in * m_coeff;
            return (delayed + feedforward);
        }

        void clearBuffers() noexcept override {
            m_delayLine.reset();
            m_lfo.setPhase(0.0f);
            m_samplesUntilUpdate = 0;
        }

    private:
        const int m_updateRate{ 100 };
        int m_samplesUntilUpdate{ 0 };
        float m_excursionSeconds{ 0.0f };
        SDSP::Oscillators::SDSPOscillator m_lfo;
        juce::SmoothedValue<float> m_smoothedDelayTimeSamples;
    };

    class [[maybe_unused]] APF2ndOrder {
    public:
        [[maybe_unused]] float tap(const int index) {
            return m_delayLine.popSample(0, static_cast<float>(index), false);
        }

        void prepareToPlay(int samplesPerBlockExpected, double sampleRate) {
            m_delayLine.prepare({sampleRate, static_cast<juce::uint32>(samplesPerBlockExpected), 1});
            m_delayLine.setMaximumDelayInSamples(static_cast<int>(sampleRate / 2.0));
            m_delayLine.setDelay(static_cast<float>(m_delayTimeSamples));
            m_hasBeenPrepared = true;
        }
        float processSample(float x) noexcept {
            auto delayed = m_delayLine.popSample(0);
            m_delayLine.pushSample(0, x + (delayed * m_coeff));
            delayed *= (1 - (m_coeff * m_coeff));
            auto feedforward = x * -m_coeff;
            auto out = delayed + feedforward;
            out *= (1 - (m_coeff * m_coeff));
            feedforward = x * -m_coeff;
            return out + feedforward;
        }

        [[maybe_unused]] void setDelayTimeSamples(int newDelayTimeSamples) {
            m_delayTimeSamples = newDelayTimeSamples;
            if(m_hasBeenPrepared) {
                m_delayLine.setDelay(static_cast<float>(m_delayTimeSamples));
            }
        }

        [[maybe_unused]] void setCoeff(float newCoeff) {
            m_coeff = newCoeff;
        }

    private:
        bool m_hasBeenPrepared{ false };
        float m_coeff{};
        int m_delayTimeSamples{};
        juce::dsp::DelayLine<float> m_delayLine;
    };


}
