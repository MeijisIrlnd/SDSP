#pragma once 
#include <juce_core/juce_core.h>
#include "../Macros.h"
namespace SDSP 
{ 
    // Mostly lifted from here: https://github.com/stekyne/PhaseVocoder/blob/master/DSP/BlockCircularBuffer.h
    // No license was provided, but please do check out his github https://github.com/stekyne
    template <typename T>
    class BlockCircularBuffer
    {
    public:
        BlockCircularBuffer() = default;
        BlockCircularBuffer(size_t size) {
            setSize(size, true);
        }

        SDSP_INLINE void setReadHopSize(size_t newSize) {
            m_readHopSize = newSize;
        }

        SDSP_INLINE void setWriteHopSize(size_t newSize) {
            m_writeHopSize = newSize;
        }

        SDSP_INLINE const size_t getReadHopSize() { return m_readHopSize; }
        SDSP_INLINE const size_t getWriteHopSize() { return m_writeHopSize; }
        SDSP_INLINE const size_t getReadIndex() { return m_readIndex; }
        SDSP_INLINE const size_t getWriteIndex() { return m_writeIndex; }

        SDSP_INLINE void setSize(size_t newSize, bool clear = false)
        {
            if (newSize == m_length) {
                if (clear) {
                    reset();
                }
                return;
            }
            m_block.allocate(newSize, clear);
            m_length = newSize;
            m_writeIndex = 0;
            m_readIndex = 0;

        }

        void reset() {
            m_block.clear(m_length);
            m_writeIndex = 0;
            m_readIndex = 0;
        }

        // read from the internal buffer into a destination buffer, wrapping around in the internal buffer if necessary
        SDSP_INLINE void read(T* dest, const size_t destLength) noexcept
        {
            const auto firstReadAmt = m_readIndex + destLength >= m_length ? m_length - m_readIndex : destLength;
            const auto internalBuffer = m_block.getData();
            std::memcpy(dest, internalBuffer + m_readIndex, sizeof(T) * firstReadAmt);
            if (firstReadAmt < destLength) {
                //ulong long = ptr I think...
                std::memcpy(dest + firstReadAmt, internalBuffer, sizeof(T) * static_cast<unsigned long long>(destLength) - firstReadAmt);
            }
            m_readIndex += m_readHopSize != 0 ? m_readHopSize : destLength;
            m_readIndex %= m_length;
        }

        SDSP_INLINE void write(const T* source, size_t sourceLength)
        {
            const auto firstWriteAmt = m_writeIndex + sourceLength >= m_length ? m_length - m_writeIndex : sourceLength;
            auto internalBuffer = m_block.getData();
            std::memcpy(internalBuffer + m_writeIndex, source, sizeof(T) * firstWriteAmt);
            if (firstWriteAmt < sourceLength) {
                std::memcpy(internalBuffer, source + firstWriteAmt, sizeof(T) * static_cast<unsigned long long>(sourceLength) - firstWriteAmt);
            }
            m_writeIndex += m_writeHopSize != 0 ? m_writeHopSize : sourceLength;
            m_writeIndex %= m_length;
            m_latestDataIndex = m_writeIndex + sourceLength % m_length;
        }

        SDSP_INLINE void overlapWrite(const T* source, size_t sourceLength)
        {
            const int writeIndexDelta = getIndexDelta(m_writeIndex, m_latestDataIndex, m_length);
            const int overlapSampleCount = sourceLength - m_writeHopSize;
            const auto overlapAmount = std::min(writeIndexDelta, overlapSampleCount);
            auto tempWriteIndex = m_writeIndex;
            auto firstWriteAmt = m_writeIndex + overlapAmount > m_length ? m_length - m_writeIndex : overlapAmount;
            auto* internalBuffer = m_block.getData();
            juce::FloatVectorOperations::add(internalBuffer + m_writeIndex, source, firstWriteAmt);
            if (firstWriteAmt < overlapAmount) {
                juce::FloatVectorOperations::add(internalBuffer, source + firstWriteAmt, overlapAmount - firstWriteAmt);
            }
            tempWriteIndex += overlapAmount;
            tempWriteIndex %= m_length;
            const auto remaining = sourceLength - overlapAmount;
            firstWriteAmt = tempWriteIndex + remaining > m_length ? m_length - tempWriteIndex : remaining;
            std::memcpy(internalBuffer + tempWriteIndex, source + overlapAmount, sizeof(T) * firstWriteAmt);
            if (firstWriteAmt < remaining) {
                std::memcpy(internalBuffer, source + overlapAmount + firstWriteAmt, sizeof(T) * (remaining - static_cast<unsigned long long>(firstWriteAmt)));
            }
            m_latestDataIndex = (m_writeIndex + sourceLength) % m_length;
            m_writeIndex += m_writeHopSize;
            m_writeIndex %= m_length;
        }

    private:
        SDSP_INLINE int getIndexDelta(int i1, int i2, int bufferLen) noexcept {
            return (i1 <= i2) ? i2 - i1 : bufferLen - i1 + i2;
        }
        size_t m_length{ 0 };
        size_t m_writeIndex{ 0 }, m_readIndex{ 0 };
        size_t m_readHopSize{ 0 }, m_writeHopSize{ 0 };
        size_t m_latestDataIndex{ 0 };
        juce::HeapBlock<T> m_block;
    };
}