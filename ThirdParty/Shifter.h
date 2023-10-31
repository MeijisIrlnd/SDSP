//
// Created by Syl on 28/03/2023.
//
#pragma once
#include <juce_core/juce_core.h>
#include "../Macros.h"
#include "SignalsmithDSP.h"
#include <array>
using SSpectrum = signalsmith::spectral::STFT<float>::Spectrum;
namespace SDSP
{
    class Shifter
    {
    public:
        struct StereoSample{
            std::array<float, 2> samples{ 0.0f, 0.0f};
            [[nodiscard]] float left() noexcept {
                return samples[0];
            }

            [[nodiscard]]float right() noexcept {
                return samples[1];
            }

            [[nodiscard]] float& operator[](size_t channel) {
                return samples[channel];
            }
        };

        Shifter(int fftSize, int hopSize) : m_fftSize(fftSize), m_hopSize(hopSize) {
            for(size_t i = 0; i < 2; ++i) {
                m_accumulator[i].resize(fftSize);
                const auto nBins = m_fftSize / 2;
                m_analysisFreqs[i].resize(nBins);
                m_analysisMags[i].resize(nBins);
                m_synthesisFreqs[i].resize(nBins);
                m_synthesisMags[i].resize(nBins);
                m_phaseAccumulator[i].resize(nBins);
                m_prevPhases[i].resize(nBins);
            }
        }

        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate) {
            m_sampleRate = sampleRate;
            m_stft.resize(2, m_fftSize, m_hopSize);
            m_stft.reset();

            for(size_t i = 0; i < 2; ++i) {
                std::fill(m_accumulator[i].begin(), m_accumulator[i].end(), 0.0f);
                std::fill(m_analysisFreqs[i].begin(), m_analysisFreqs[i].end(), 0.0f);
                std::fill(m_analysisMags[i].begin(), m_analysisMags[i].end(), 0.0f);
                std::fill(m_synthesisFreqs[i].begin(), m_synthesisFreqs[i].end(), 0.0f);
                std::fill(m_synthesisMags[i].begin(), m_synthesisMags[i].end(), 0.0f);
                std::fill(m_phaseAccumulator[i].begin(), m_phaseAccumulator[i].end(), 0.0f);
                std::fill(m_prevPhases[i].begin(), m_prevPhases[i].end(), 0.0f);
            }

        }

        StereoSample processSample(StereoSample& sample) noexcept {
            for(size_t i = 0; i < 2; ++i) {
                pushToAccumulator(i, sample[i]);
            }
            m_stft.ensureValid([&](int) {
                m_stft.analyse(m_accumulator);
                stftCallback(m_stft.spectrum);
            });
            auto res = m_stft.at(0);
            ++m_stft;
            return {static_cast<float>(res[0]), static_cast<float>(res[1])};
        }

        SDSP_INLINE void setShift(float newShift) noexcept {
            m_shift = newShift;
        }

    private:
        struct MagPhasePair {
            float mag{}, phase{};
        };

        void pushToAccumulator(size_t channel, float toPush) {
            m_accumulator[channel].emplace_back(toPush);
            m_accumulator[channel].erase(m_accumulator[channel].begin(), m_accumulator[channel].begin() + 1);
        }

        void stftCallback(SSpectrum& spectrum) {
            MagPhasePair mpPair{0.0f, 0.0f};
            const auto freqPerBin = static_cast<float>(m_sampleRate) / static_cast<float>(m_fftSize);
            const auto nBins = m_stft.bands();
            const auto oversamplingFactor = static_cast<float>(m_fftSize) / static_cast<float>(m_hopSize);
            const auto expectedPhaseDelta = juce::MathConstants<float>::twoPi * (static_cast<float>(m_hopSize) / static_cast<float>(m_fftSize));
            for(auto channel = 0; channel < 2; ++channel) {
                auto* data = spectrum[channel];
                for(auto bin = 0; bin < nBins; ++bin) {
                    auto currentBin = data[bin];
                    // SMB shift multiplies this by 2..
                    mpPair.mag = std::sqrt(currentBin.real() * currentBin.real() + currentBin.imag() * currentBin.imag());
                    mpPair.phase = std::atan2(currentBin.imag(), currentBin.real());

                    auto temp = mpPair.phase - m_prevPhases[static_cast<size_t>(channel)][static_cast<size_t>(bin)];
                    m_prevPhases[static_cast<size_t>(channel)][static_cast<size_t>(bin)] = mpPair.phase;
                    temp -= static_cast<float>(bin) * expectedPhaseDelta;

                    // map delta phase to -pi to pi

                    auto qpd = temp / juce::MathConstants<float>::pi;
                    if (qpd >= 0) qpd += static_cast<int>(qpd) & 1;
                    else qpd -= static_cast<int>(qpd) & 1;
                    temp -= juce::MathConstants<float>::pi * (double) qpd;

                    // deviation from bin freq from interval..
                    temp = oversamplingFactor * temp / juce::MathConstants<float>::twoPi;
                    // compute partial's true freq..
                    temp = static_cast<float>(bin) * freqPerBin + temp * freqPerBin;
                    m_analysisMags[channel][bin] = mpPair.mag;
                    m_analysisFreqs[channel][bin] = mpPair.phase;
                }

                // Shift ..
                for(size_t bin = 0; bin < static_cast<size_t>(nBins); ++bin) {
                    auto index = static_cast<size_t>(static_cast<float>(bin) * m_shift);
                    if(index < nBins) {
                        m_synthesisMags[static_cast<size_t>(channel)][index] = m_analysisMags[static_cast<size_t>(channel)][bin];
                        m_synthesisFreqs[static_cast<size_t>(channel)][index] = m_analysisFreqs[static_cast<size_t>(channel)][bin] * m_shift;
                    }
                }

                // Synthesis now...
                for(auto bin = 0; bin < nBins; ++bin) {
                    auto mag = m_synthesisMags[static_cast<size_t>(channel)][static_cast<size_t>(bin)];
                    auto temp = m_synthesisFreqs[static_cast<size_t>(channel)][static_cast<size_t>(bin)];

                    temp -= static_cast<float>(bin) * freqPerBin;
                    temp /= freqPerBin;
                    temp = juce::MathConstants<float>::twoPi * (temp / oversamplingFactor);
                    temp += static_cast<float>(bin) * expectedPhaseDelta;

                    m_phaseAccumulator[static_cast<size_t>(channel)][static_cast<size_t>(bin)] += temp;
                    auto phase = m_phaseAccumulator[static_cast<size_t>(channel)][static_cast<size_t>(bin)];
                    spectrum[channel][static_cast<size_t>(bin)] = {mag * std::cos(phase), mag * std::sin(phase)};
                }
            }
        }

        double m_sampleRate{ 44100.0 };
        signalsmith::spectral::STFT<float> m_stft;
        std::array<std::vector<float>, 2> m_accumulator;
        std::array<std::vector<float>, 2> m_prevPhases;
        std::array<std::vector<float>, 2> m_analysisFreqs, m_analysisMags, m_synthesisMags, m_synthesisFreqs, m_phaseAccumulator;
        int m_fftSize{}, m_hopSize{} ;

        float m_shift{ 1.0f };
    };
}