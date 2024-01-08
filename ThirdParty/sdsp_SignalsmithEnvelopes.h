#ifndef SDSP_SIGNALSMITHENVELOPES_H
#define SDSP_SIGNALSMITHENVELOPES_H
#include <vector>
namespace SDSP::Signalsmith::Envelopes {
    /** Variable-width rectangular sum */
    template <typename Sample = double>
    class BoxSum {
        int bufferLength, index;
        std::vector<Sample> buffer;
        Sample sum = 0, wrapJump = 0;

    public:
        BoxSum(int maxLength) {
            resize(maxLength);
        }

        /// Sets the maximum size (and reset contents)
        void resize(int maxLength) {
            bufferLength = maxLength + 1;
            buffer.resize(bufferLength);
            if (maxLength != 0) buffer.shrink_to_fit();
            reset();
        }

        /// Resets (with an optional "fill" value)
        void reset(Sample value = Sample()) {
            index = 0;
            sum = 0;
            for (size_t i = 0; i < buffer.size(); ++i) {
                buffer[i] = sum;
                sum += value;
            }
            wrapJump = sum;
            sum = 0;
        }

        Sample read(int width) {
            int readIndex = index - width;
            double result = sum;
            if (readIndex < 0) {
                result += wrapJump;
                readIndex += bufferLength;
            }
            return static_cast<Sample>(result) - buffer[readIndex];
        }

        void write(Sample value) {
            ++index;
            if (index == bufferLength) {
                index = 0;
                wrapJump = sum;
                sum = 0;
            }
            sum += value;
            buffer[index] = sum;
        }

        Sample readWrite(Sample value, int width) {
            write(value);
            return read(width);
        }
    };

    /** Rectangular moving average filter (FIR).
            \diagram{box-filter-example.svg}
            A filter of length 1 has order 0 (i.e. does nothing). */
    template <typename Sample = double>
    class BoxFilter {
        BoxSum<Sample> boxSum;
        int _length, _maxLength;
        Sample multiplier;

    public:
        BoxFilter(int maxLength) : boxSum(maxLength) {
            resize(maxLength);
        }
        /// Sets the maximum size (and current size, and resets)
        void resize(int maxLength) {
            _maxLength = maxLength;
            boxSum.resize(maxLength);
            set(maxLength);
        }
        /// Sets the current size (expanding/allocating only if needed)
        void set(int length) {
            _length = length;
            multiplier = Sample(1) / length;
            if (length > _maxLength) resize(length);
        }

        /// Resets (with an optional "fill" value)
        void reset(Sample fill = Sample()) {
            boxSum.reset(fill);
        }

        Sample operator()(Sample v) {
            return boxSum.readWrite(v, _length) * multiplier;
        }
    };

    template <typename Sample>
    class PeakHold {
        static constexpr Sample lowest = std::numeric_limits<Sample>::lowest();
        int bufferMask;
        std::vector<Sample> buffer;
        int backIndex = 0, middleStart = 0, workingIndex = 0, middleEnd = 0, frontIndex = 0;
        Sample frontMax = lowest, workingMax = lowest, middleMax = lowest;

    public:
        PeakHold(int maxLength) {
            resize(maxLength);
        }
        int size() {
            return frontIndex - backIndex;
        }
        void resize(int maxLength) {
            int bufferLength = 1;
            while (bufferLength < maxLength) bufferLength *= 2;
            buffer.resize(bufferLength);
            bufferMask = bufferLength - 1;

            frontIndex = backIndex + maxLength;
            reset();
        }
        void reset(Sample fill = lowest) {
            int prevSize = size();
            buffer.assign(buffer.size(), fill);
            frontMax = workingMax = middleMax = lowest;
            middleEnd = workingIndex = frontIndex = 0;
            middleStart = middleEnd - (prevSize / 2);
            backIndex = frontIndex - prevSize;
        }
        /** Sets the size immediately.
        Must be `0 <= newSize <= maxLength` (see constructor and `.resize()`).

        Shrinking doesn't destroy information, and if you expand again (with `preserveCurrentPeak=false`), you will get the same output as before shrinking.  Expanding when `preserveCurrentPeak` is enabled is destructive, re-writing its history such that the current output value is unchanged.*/
        void set(int newSize, bool preserveCurrentPeak = false) {
            while (size() < newSize) {
                Sample& backPrev = buffer[backIndex & bufferMask];
                --backIndex;
                Sample& back = buffer[backIndex & bufferMask];
                back = preserveCurrentPeak ? backPrev : std::max(back, backPrev);
            }
            while (size() > newSize) {
                pop();
            }
        }

        void push(Sample v) {
            buffer[frontIndex & bufferMask] = v;
            ++frontIndex;
            frontMax = std::max(frontMax, v);
        }
        void pop() {
            if (backIndex == middleStart) {
                // Move along the maximums
                workingMax = lowest;
                middleMax = frontMax;
                frontMax = lowest;

                int prevFrontLength = frontIndex - middleEnd;
                int prevMiddleLength = middleEnd - middleStart;
                if (prevFrontLength <= prevMiddleLength + 1) {
                    // Swap over simply
                    middleStart = middleEnd;
                    middleEnd = frontIndex;
                    workingIndex = middleEnd;
                } else {
                    // The front is longer than the middle - only happens if unbalanced
                    // We don't move *all* of the front over, keeping half the surplus in the front
                    int middleLength = (frontIndex - middleStart) / 2;
                    middleStart = middleEnd;
                    middleEnd += middleLength;

                    // Working index is close enough that it will be finished by the time the back is empty
                    int backLength = middleStart - backIndex;
                    int workingLength = std::min(backLength, middleEnd - middleStart);
                    workingIndex = middleStart + workingLength;

                    // Since the front was not completely consumed, we re-calculate the front's maximum
                    for (int i = middleEnd; i != frontIndex; ++i) {
                        frontMax = std::max(frontMax, buffer[i & bufferMask]);
                    }
                    // The index might not start at the end of the working block - compute the last bit immediately
                    for (int i = middleEnd - 1; i != workingIndex - 1; --i) {
                        buffer[i & bufferMask] = workingMax = std::max(workingMax, buffer[i & bufferMask]);
                    }
                }

                // Is the new back (previous middle) empty? Only happens if unbalanced
                if (backIndex == middleStart) {
                    // swap over again (front's empty, no change)
                    workingMax = lowest;
                    middleMax = frontMax;
                    frontMax = lowest;
                    middleStart = workingIndex = middleEnd;

                    if (backIndex == middleStart) {
                        --backIndex; // Only happens if you pop from an empty list - fail nicely
                    }
                }

                buffer[frontIndex & bufferMask] = lowest; // In case of length 0, when everything points at this value
            }

            ++backIndex;
            if (workingIndex != middleStart) {
                --workingIndex;
                buffer[workingIndex & bufferMask] = workingMax = std::max(workingMax, buffer[workingIndex & bufferMask]);
            }
        }
        Sample read() {
            Sample backMax = buffer[backIndex & bufferMask];
            return std::max(backMax, std::max(middleMax, frontMax));
        }

        // For simple use as a constant-length filter
        Sample operator()(Sample v) {
            push(v);
            pop();
            return read();
        }
    };
} // namespace SDSP::Signalsmith::Envelopes
#endif