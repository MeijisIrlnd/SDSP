//
// Created by Syl on 22/03/2023.
//

/****************************************************************************
*
* NAME: smbPitchShift.cpp
* VERSION: 1.2
* HOME URL: http://blogs.zynaptiq.com/bernsee
* KNOWN BUGS: none
*
* SYNOPSIS: Routine for doing pitch shifting while maintaining
* duration ustd::sing the Short Time Fourier Transform.
*
* DESCRIPTION: The routine takes a pitchShift factor value which is between 0.5
* (one octave down) and 2. (one octave up). A value of exactly 1 does not change
* the pitch. numSampsToProcess tells the routine how many samples in indata[0...
* numSampsToProcess-1] should be pitch shifted and moved to outdata[0 ...
* numSampsToProcess-1]. The two buffers can be identical (ie. it can process the
* data in-place). fftFrameSize defines the FFT frame size used for the
* processtd::sing. Typical values are 1024, 2048 and 4096. It may be any value <=
* MAX_FRAME_LENGTH but it MUST be a power of 2. osamp is the STFT
* oversampling factor which also determines the overlap between adjacent STFT
* frames. It should at least be 4 for moderate scaling ratios. A value of 32 is
* recommended for best quality. sampleRate takes the sample rate for the signal
* in unit Hz, ie. 44100 for 44.1 kHz audio. The data passed to the routine in
* indata[] should be in the range [-1.0, 1.0), which is also the output range
* for the data, make sure you scale the data accordingly (for 16bit signed integers
* you would have to divide (and multiply) by 32768).
*
* COPYRIGHT 1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
*
* 						The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies.
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wol.htm for more information.
*
*****************************************************************************/
#pragma once
#include "SignalsmithDSP.h"
#include <juce_dsp/juce_dsp.h>
#include <string>>
#include <cmath>
#include <array>
#include "../Macros.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define MAX_FRAME_LENGTH 8192
namespace SMB {
    static void smbFft(float *fftBuffer, long fftFrameSize, long sign);

    static double smbAtan2(double x, double y);


// -----------------------------------------------------------------------------------------------------------------
    class SMBPitchShift
    {
    private:
        float m_shift{ 1.0f };
        double m_sampleRate{ 44100.0 };
        std::array<float, MAX_FRAME_LENGTH> m_inFifo{}, m_outFifo{}, m_analysisFrequencies{}, m_analysisMagnitudes{}, m_synthesisFrequencies{}, m_synthesisMagnitudes{};
        std::array<float, 2 * MAX_FRAME_LENGTH> m_fftWorkspace{}, m_outputAccumulator{};
        std::array<float, MAX_FRAME_LENGTH / 2 + 1> m_lastPhase{}, m_sumPhase{};
        std::vector<std::complex<float> > m_spectrum{};
        long m_rOver{ false };
        long m_fftFrameSize{}, m_oversamplingFactor{};
        signalsmith::fft::RealFFT<float> m_fft;
    public:
        explicit SMBPitchShift(long fftFrameSize, long oversamplingFactor) : m_fftFrameSize(fftFrameSize), m_oversamplingFactor(oversamplingFactor) {
            std::memset(m_inFifo.data(), 0, MAX_FRAME_LENGTH * sizeof(float));
            std::memset(m_outFifo.data(), 0, MAX_FRAME_LENGTH * sizeof(float));
            std::memset(m_fftWorkspace.data(), 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
            std::memset(m_lastPhase.data(), 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
            std::memset(m_sumPhase.data(), 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
            std::memset(m_outputAccumulator.data(), 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
            std::memset(m_analysisFrequencies.data(), 0, MAX_FRAME_LENGTH * sizeof(float));
            std::memset(m_analysisMagnitudes.data(), 0, MAX_FRAME_LENGTH * sizeof(float));
            m_spectrum.resize(static_cast<size_t>((fftFrameSize) / 2) + 1);
            std::fill(m_spectrum.begin(), m_spectrum.end(), std::complex<float>{0.0f, 0.0f});
            m_fft.setSize(static_cast<size_t>(fftFrameSize));
        }

        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate) {
            m_sampleRate = sampleRate;
        }

        float processSample(float x) noexcept {
            auto out{ 0.0f };
            double magn, phase, tmp, window, real, imag;
            double freqPerBin, expct;
            long i, k, qpd, index, inFifoLatency, stepSize, fftFrameSize2;

            /* set up some handy variables */
            fftFrameSize2 = m_fftFrameSize / 2;
            stepSize = m_fftFrameSize / m_oversamplingFactor;
            freqPerBin = m_sampleRate / (double) m_fftFrameSize;
            expct = 2. * M_PI * (double) stepSize / (double) m_fftFrameSize;
            inFifoLatency = m_fftFrameSize - stepSize;
            if (m_rOver == false) m_rOver = inFifoLatency;


            /* As long as we have not yet collected enough data just read in */
            m_inFifo[m_rOver] = x;
            out = m_outFifo[m_rOver - inFifoLatency];
            m_rOver++;

            /* now we have enough data for processtd::sing */
            if (m_rOver >= m_fftFrameSize) {
                m_rOver = inFifoLatency;

                /* do windowing and re,im interleave */
                for (k = 0; k < m_fftFrameSize; k++) {
                    window = -.5 * juce::dsp::FastMathApproximations::cos(2. * M_PI * (double) k / (double) m_fftFrameSize) + .5;
                    m_fftWorkspace[2 * k] = m_inFifo[k] * window;
                    m_fftWorkspace[2 * k + 1] = 0.;
                }


                /* ***************** ANALYSIS ******************* */
                /* do transform */
                // replace this with Signalsmith's fft?

                smbFft(m_fftWorkspace.data(), m_fftFrameSize, -1);
                /* this is the analysis step */
                for (k = 0; k <= fftFrameSize2; k++) {

                    /* de-interlace FFT buffer */
                    real = m_fftWorkspace[2 * k];
                    imag = m_fftWorkspace[2 * k + 1];

                    /* compute magnitude and phase */
                    magn = 2. * sqrt(real * real + imag * imag);
                    phase = atan2(imag, real);

                    /* compute phase difference */
                    tmp = phase - m_lastPhase[k];
                    m_lastPhase[k] = phase;

                    /* subtract expected phase difference */
                    tmp -= (double) k * expct;

                    /* map delta phase into +/- Pi interval */
                    qpd = tmp / M_PI;
                    if (qpd >= 0) qpd += qpd & 1;
                    else qpd -= qpd & 1;
                    tmp -= M_PI * (double) qpd;

                    /* get deviation from bin frequency from the +/- Pi interval */
                    tmp = m_oversamplingFactor * tmp / (2. * M_PI);

                    /* compute the k-th partials' true frequency */
                    tmp = (double) k * freqPerBin + tmp * freqPerBin;

                    /* store magnitude and true frequency in analysis arrays */
                    m_analysisMagnitudes[k] = magn;
                    m_analysisFrequencies[k] = tmp;

                }

                /* ***************** PROCESstd::sinG ******************* */
                /* this does the actual pitch shifting */
                memset(m_synthesisMagnitudes.data(), 0, m_fftFrameSize * sizeof(float));
                memset(m_synthesisFrequencies.data(), 0, m_fftFrameSize * sizeof(float));
                for (k = 0; k <= fftFrameSize2; k++) {
                    index = k * m_shift;
                    if (index <= fftFrameSize2) {
                        m_synthesisMagnitudes[index] += m_analysisMagnitudes[k];
                        m_synthesisFrequencies[index] = m_analysisFrequencies[k] * m_shift;
                    }
                }

                /* ***************** SYNTHESIS ******************* */
                /* this is the synthesis step */
                for (k = 0; k <= fftFrameSize2; k++) {

                    /* get magnitude and true frequency from synthesis arrays */
                    magn = m_synthesisMagnitudes[k];
                    tmp = m_synthesisFrequencies[k];

                    /* subtract bin mid frequency */
                    tmp -= (double) k * freqPerBin;

                    /* get bin deviation from freq deviation */
                    tmp /= freqPerBin;

                    /* take osamp into account */
                    tmp = 2. * M_PI * tmp / m_oversamplingFactor;

                    /* add the overlap phase advance back in */
                    tmp += (double) k * expct;

                    /* accumulate delta phase to get bin phase */
                    m_sumPhase[k] += tmp;
                    phase = m_sumPhase[k];

                    /* get real and imag part and re-interleave */
                    m_fftWorkspace[2 * k] = magn * juce::dsp::FastMathApproximations::cos(phase);
                    m_fftWorkspace[2 * k + 1] = magn * std::sin(phase);
                }

                /* zero negative frequencies */
                for (k = m_fftFrameSize + 2; k < 2 * m_fftFrameSize; k++) m_fftWorkspace[k] = 0.;

                /* do inverse transform */
                smbFft(m_fftWorkspace.data(), m_fftFrameSize, 1);

                /* do windowing and add to output accumulator */
                for (k = 0; k < m_fftFrameSize; k++) {
                    window = -.5 * juce::dsp::FastMathApproximations::cos(2. * M_PI * (double) k / (double) m_fftFrameSize) + .5;
                    m_outputAccumulator[k] += 2. * window * m_fftWorkspace[2 * k] / (fftFrameSize2 * m_oversamplingFactor);
                }
                for (k = 0; k < stepSize; k++) m_outFifo[k] = m_outputAccumulator[k];

                /* shift accumulator */
                memmove(m_outputAccumulator.data(), m_outputAccumulator.data() + stepSize,
                        m_fftFrameSize * sizeof(float));

                /* move input FIFO */
                for (k = 0; k < inFifoLatency; k++) m_inFifo[k] = m_inFifo[k + stepSize];
            }
            return out;

        }

        float processSampleSSFFT(float x) noexcept {
            auto out{ 0.0f };
            double magn, phase, tmp, window, real, imag;
            double freqPerBin, expct;
            long i, k, qpd, index, inFifoLatency, stepSize, fftFrameSize2;

            /* set up some handy variables */
            fftFrameSize2 = m_fftFrameSize / 2;
            stepSize = m_fftFrameSize / m_oversamplingFactor;
            freqPerBin = m_sampleRate / (double) m_fftFrameSize;
            expct = 2. * M_PI * (double) stepSize / (double) m_fftFrameSize;
            inFifoLatency = m_fftFrameSize - stepSize;
            if (m_rOver == false) m_rOver = inFifoLatency;


            /* As long as we have not yet collected enough data just read in */
            m_inFifo[m_rOver] = x;
            out = m_outFifo[m_rOver - inFifoLatency];
            m_rOver++;

            /* now we have enough data for processtd::sing */
            if (m_rOver >= m_fftFrameSize) {
                m_rOver = inFifoLatency;

                // could remove this and use Signalsmith's windowed fft class instead..
                /* do windowing and re,im interleave */
                for (k = 0; k < m_fftFrameSize; k++) {
                    window = -.5 * juce::dsp::FastMathApproximations::cos(2. * M_PI * (double) k / (double) m_fftFrameSize) + .5;
                    m_fftWorkspace[2 * k] = m_inFifo[k] * window;
                    m_fftWorkspace[2 * k + 1] = 0.;
                }


                /* ***************** ANALYSIS ******************* */
                /* do transform */
                // replace this with Signalsmith's fft?
                m_fft.fft(m_fftWorkspace, m_spectrum);
                /* this is the analysis step */
                for (k = 0; k <= fftFrameSize2; k++) {
                    /* compute magnitude and phase */
                    auto currentBin = m_spectrum[k];
                    magn = 2 * std::sqrt(std::pow(currentBin.real(), 2) + std::pow(currentBin.imag(), 2));
                    phase = atan2(currentBin.imag(), currentBin.real());

                    /* compute phase difference */
                    tmp = phase - m_lastPhase[k];
                    m_lastPhase[k] = phase;

                    /* subtract expected phase difference */
                    tmp -= (double) k * expct;

                    /* map delta phase into +/- Pi interval */
                    qpd = tmp / M_PI;
                    if (qpd >= 0) qpd += qpd & 1;
                    else qpd -= qpd & 1;
                    tmp -= M_PI * (double) qpd;

                    /* get deviation from bin frequency from the +/- Pi interval */
                    tmp = m_oversamplingFactor * tmp / (2. * M_PI);

                    /* compute the k-th partials' true frequency */
                    tmp = (double) k * freqPerBin + tmp * freqPerBin;

                    /* store magnitude and true frequency in analysis arrays */
                    m_analysisMagnitudes[k] = magn;
                    m_analysisFrequencies[k] = tmp;

                }

                /* ***************** PROCESstd::sinG ******************* */
                /* this does the actual pitch shifting */
                memset(m_synthesisMagnitudes.data(), 0, m_fftFrameSize * sizeof(float));
                memset(m_synthesisFrequencies.data(), 0, m_fftFrameSize * sizeof(float));
                for (k = 0; k <= fftFrameSize2; k++) {
                    index = k * m_shift;
                    if (index <= fftFrameSize2) {
                        m_synthesisMagnitudes[index] += m_analysisMagnitudes[k];
                        m_synthesisFrequencies[index] = m_analysisFrequencies[k] * m_shift;
                    }
                }

                /* ***************** SYNTHESIS ******************* */
                /* this is the synthesis step */
                for (k = 0; k <= fftFrameSize2; k++) {

                    /* get magnitude and true frequency from synthesis arrays */
                    magn = m_synthesisMagnitudes[k];
                    tmp = m_synthesisFrequencies[k];

                    /* subtract bin mid frequency */
                    tmp -= (double) k * freqPerBin;

                    /* get bin deviation from freq deviation */
                    tmp /= freqPerBin;

                    /* take osamp into account */
                    tmp = 2. * M_PI * tmp / m_oversamplingFactor;

                    /* add the overlap phase advance back in */
                    tmp += (double) k * expct;

                    /* accumulate delta phase to get bin phase */
                    m_sumPhase[k] += tmp;
                    phase = m_sumPhase[k];

                    /* get real and imag part and re-interleave */
                    // don't fucking reinterleave lmaooo
                    m_spectrum[k] = { static_cast<float>(magn * juce::dsp::FastMathApproximations::cos(phase)), static_cast<float>(magn * std::sin(phase)) };
                    //m_fftWorkspace[2 * k] = magn * juce::dsp::FastMathApproximations::cos(phase);
                    //m_fftWorkspace[2 * k + 1] = magn * std::sin(phase);
                }

                /* zero negative frequencies */
                for (k = m_fftFrameSize + 2; k < 2 * m_fftFrameSize; k++) m_fftWorkspace[k] = 0.;

                m_fft.ifft(m_spectrum, m_fftWorkspace);

                /* do windowing and add to output accumulator */
                for (k = 0; k < m_fftFrameSize; k++) {
                    window = -.5 * juce::dsp::FastMathApproximations::cos(2. * M_PI * (double) k / (double) m_fftFrameSize) + .5;
                    m_outputAccumulator[k] += 2. * window * m_fftWorkspace[2 * k] / (fftFrameSize2 * m_oversamplingFactor);
                }
                for (k = 0; k < stepSize; k++) m_outFifo[k] = m_outputAccumulator[k];

                /* shift accumulator */
                memmove(m_outputAccumulator.data(), m_outputAccumulator.data() + stepSize,
                        m_fftFrameSize * sizeof(float));

                /* move input FIFO */
                for (k = 0; k < inFifoLatency; k++) m_inFifo[k] = m_inFifo[k + stepSize];
            }
            return out;

        }
        SDSP_INLINE void setShift(float newShift) noexcept {
            m_shift = newShift;
        }
    };



// -----------------------------------------------------------------------------------------------------------------


    static inline void smbFft(float *fftBuffer, long fftFrameSize, long sign)
/*
	FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
	Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
	time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
	and returns the juce::dsp::FastMathApproximations::costd::sine and std::sine parts in an interleaved manner, ie.
	fftBuffer[0] = juce::dsp::FastMathApproximations::cosPart[0], fftBuffer[1] = std::sinPart[0], asf. fftFrameSize
	must be a power of 2. It expects a complex input signal (see footnote 2),
	ie. when working with 'common' audio signals our input signal has to be
	passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
	of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
    {
        float wr, wi, arg, *p1, *p2, temp;
        float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
        long i, bitm, j, le, le2, k;

        for (i = 2; i < 2 * fftFrameSize - 2; i += 2) {
            for (bitm = 2, j = 0; bitm < 2 * fftFrameSize; bitm <<= 1) {
                if (i & bitm) j++;
                j <<= 1;
            }
            if (i < j) {
                p1 = fftBuffer + i;
                p2 = fftBuffer + j;
                temp = *p1;
                *(p1++) = *p2;
                *(p2++) = temp;
                temp = *p1;
                *p1 = *p2;
                *p2 = temp;
            }
        }
        for (k = 0, le = 2; k < (long) (log(fftFrameSize) / log(2.) + .5); k++) {
            le <<= 1;
            le2 = le >> 1;
            ur = 1.0;
            ui = 0.0;
            arg = M_PI / (le2 >> 1);
            wr = juce::dsp::FastMathApproximations::cos(arg);
            wi = sign * std::sin(arg);
            for (j = 0; j < le2; j += 2) {
                p1r = fftBuffer + j;
                p1i = p1r + 1;
                p2r = p1r + le2;
                p2i = p2r + 1;
                for (i = j; i < 2 * fftFrameSize; i += le) {
                    tr = *p2r * ur - *p2i * ui;
                    ti = *p2r * ui + *p2i * ur;
                    *p2r = *p1r - tr;
                    *p2i = *p1i - ti;
                    *p1r += tr;
                    *p1i += ti;
                    p1r += le;
                    p1i += le;
                    p2r += le;
                    p2i += le;
                }
                tr = ur * wr - ui * wi;
                ui = ur * wi + ui * wr;
                ur = tr;
            }
        }
    }


// -----------------------------------------------------------------------------------------------------------------

/*

    12/12/02, smb

    PLEASE NOTE:

    There have been some reports on domain errors when the atan2() function was used
    as in the above code. Usually, a domain error should not interrupt the program flow
    (maybe except in Debug mode) but rather be handled "silently" and a global variable
    should be set according to this error. However, on some occasions people ran into
    this kind of scenario, so a replacement atan2() function is provided here.

    If you are experiencing domain errors and your program stops, simply replace all
    instances of atan2() with calls to the smbAtan2() function below.

*/


    static inline double smbAtan2(double x, double y) {
        double signx;
        if (x > 0.) signx = 1.;
        else signx = -1.;

        if (x == 0.) return 0.;
        if (y == 0.) return signx * M_PI / 2.;

        return atan2(x, y);
    }


// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
}// SMB