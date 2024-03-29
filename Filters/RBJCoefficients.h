/*
  ==============================================================================

    RBJCoefficients.h
    Created: 2 Jan 2022 3:35:54pm
    Author:  Syl

  ==============================================================================
*/

#pragma once
#include <juce_core/juce_core.h>
#include "../KMath.h"
#include "../Macros.h"
namespace SDSP::RBJ
    {
        SDSP_UNUSED static inline void lowpass(double* target, double sampleRate, double cutoff, double q)
        {
            const double omega = juce::MathConstants<double>::twoPi * (cutoff / sampleRate);
            const double cosOmega = std::cos(omega);
            const double alpha = std::sin(omega) / static_cast<double>(2 * q);
            target[0] = (1 - cosOmega) / 2.0; // a0
            target[1] = 1 - cosOmega; // a1
            target[2] = target[0]; // a2
            target[3] = 1 + alpha; // b0
            target[4] = -2 * cosOmega; // b1
            target[5] = 1 - alpha; // b2;
        }

        SDSP_UNUSED static inline void highpass(double* target, double sampleRate, double cutoff, double q)
        {
            const double omega = juce::MathConstants<double>::twoPi * (cutoff / sampleRate);
            const double cosOmega = std::cos(omega);
            const double alpha = std::sin(omega) / static_cast<double>(2 * q);
            target[0] = (1 + cosOmega) / 2.0; // a0
            target[1] = -1 * (1 + cosOmega); // a1
            target[2] = target[0]; // a2
            target[3] = 1 + alpha; // b0
            target[4] = -2 * cosOmega; // b1
            target[5] = 1 - alpha; // b2
        }

        SDSP_UNUSED static inline void lowShelf(double* target, double sampleRate, double centreFreq, double dbGain, double slope)
        {
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double cosOmega = std::cos(omega);
            const double sinOmega = std::sin(omega);
            const double A = std::pow(10, dbGain / 40.0);
            const double twoRootAa = sinOmega * std::pow((std::pow(A, 2) + 1) * ((1 / slope) - 1) + 2 * A, 0.5);
            target[0] = A * ((A + 1) - (A - 1) * cosOmega + twoRootAa); // a0
            target[1] = (2 * A) * ((A - 1) - (A + 1) * cosOmega); // a1
            target[2] = A * ((A + 1) - (A - 1) * cosOmega - twoRootAa); // a2
            target[3] = (A + 1) + (A - 1) * cosOmega + twoRootAa; // b0
            target[4] = -2 * ((A - 1) + (A + 1) * cosOmega); // b1
            target[5] = (A + 1) + (A - 1) * cosOmega - twoRootAa; // b2
        }

       SDSP_UNUSED static inline void highShelf(double* target, double sampleRate, double centreFreq, double dbGain, double slope)
        {
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double cosOmega = std::cos(omega);
            const double sinOmega = std::sin(omega);
            const double A = std::pow(10, dbGain / 40.0);
            const double twoRootAa = sinOmega * std::pow((std::pow(A, 2) + 1) * ((1 / slope) - 1) + 2 * A, 0.5);
            target[0] = A * ((A + 1) + (A - 1) * cosOmega + twoRootAa); // a0
            target[1] = (-2 * A) * ((A - 1) + (A + 1) * cosOmega); // a1
            target[2] = A * ((A + 1) + (A - 1) * cosOmega - twoRootAa); // a2
            target[3] = (A + 1) - (A - 1) * cosOmega + twoRootAa; // b0
            target[4] = 2 * ((A - 1) - (A + 1) * cosOmega); // b1
            target[5] = (A + 1) - (A - 1) * cosOmega - twoRootAa; // b2
        }


        // With peak gain
        SDSP_UNUSED static inline void bandpass(double* target, double sampleRate, double centreFreq, double bandwidth, double peakGain)
        {
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double sinOmega = std::sin(omega);
            const double cosOmega = std::cos(omega);
            const double ln2 = std::log(2);
            const double alpha = sinOmega * std::sinh((ln2 / 2.0) * bandwidth * (omega / sinOmega));
            target[0] = peakGain * alpha;
            target[1] = 0; // a1
            target[2] = (-1 * peakGain) * alpha; // a2
            target[3] = 1 + alpha; // b0
            target[4] = -2 * cosOmega; // b1
            target[5] = 1 - alpha; // b2
        }

        // 0db constant peak
        SDSP_UNUSED static inline void bandpass(double* target, double sampleRate, double centreFreq, double bandwidth)
        {
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double sinOmega = std::sin(omega);
            const double cosOmega = std::cos(omega);
            const double ln2 = std::log(2);
            const double alpha = sinOmega * std::sinh((ln2 / 2.0) * bandwidth * (omega / sinOmega));
            target[0] = alpha;
            target[1] = 0;
            target[2] = -1 * alpha;
            target[3] = 1 + alpha;
            target[4] = -2 * cosOmega;
            target[5] = 1 - alpha;
        }

       SDSP_UNUSED static inline void bandreject(double* target, double sampleRate, double centreFreq, double bandwidth)
        {
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double sinOmega = std::sin(omega);
            const double cosOmega = std::cos(omega);
            const double ln2 = std::log(2);
            const double alpha = sinOmega * std::sinh((ln2 / 2.0) * bandwidth * (omega / sinOmega));
            target[0] = 1;
            target[1] = -2 * cosOmega;
            target[2] = 1;
            target[3] = 1 + alpha;
            target[4] = -2 * cosOmega;
            target[5] = 1 - alpha;
        }

        // BW = width between dbGain / 2 gain frequencies
        SDSP_UNUSED static inline void bell(double* target, double sampleRate, double gainDB, double centreFreq, double bandwidth)
        {
            const double A = std::pow(10, gainDB / 40.0);
            const double omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const double sinOmega = std::sin(omega);
            const double cosOmega = std::cos(omega);
            const double ln2 = std::log(2);
            const double alpha = sinOmega * std::sinh((ln2 / 2.0) * bandwidth * (omega / sinOmega));
            target[0] = 1 + (alpha * A);
            target[1] = -2 * cosOmega;
            target[2] = 1 - (alpha * A);
            target[3] = 1 + (alpha / A);
            target[4] = -2 * cosOmega;
            target[5] = 1 - (alpha / A);
        }

        SDSP_UNUSED static inline void bellWithQ(double* target, double sampleRate, double gainDb, double centreFreq, double q) {
            const auto A = std::pow(10, gainDb / 40.0);
            const auto omega = juce::MathConstants<double>::twoPi * (centreFreq / sampleRate);
            const auto sinOmega = std::sin(omega);
            const auto cosOmega = std::cos(omega);
            const auto alpha = sinOmega / (2.0 * q);
            const auto alphaA = alpha * A;
            const auto alphaOverA = alpha / A;
            target[0] = 1 + alphaA;
            target[1] = -2 * cosOmega;
            target[2] = 1 - alphaA;
            target[3] = 1 + alphaOverA;
            target[4] = target[1];
            target[5] = 1 - alphaOverA;
        }

        SDSP_UNUSED static inline void allpass(double* target, double sampleRate, double centreFrequency, double Q)
        {
            const auto omega = juce::MathConstants<double>::twoPi * (centreFrequency / sampleRate);
            const auto alpha = std::sin(omega) / (2 * Q);
            const auto cosOmega = std::cos(omega);
            target[0] = 1 - alpha;
            target[1] = -2 * cosOmega;
            target[2] = 1 + alpha;
            target[3] = 1 + alpha;
            target[4] = -2 * cosOmega;
            target[5] = 1 - alpha;
        }
    }
