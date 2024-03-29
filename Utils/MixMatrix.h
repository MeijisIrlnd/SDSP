#pragma once
#include <cmath>

namespace SDSP {
    // Use like `Householder<double, 8>::inPlace(data)` - size must be ≥ 1
    template <typename Sample, int size>
    class [[maybe_unused]] Householder {
        static constexpr Sample multiplier{ -2.0 / size };

    public:
        static void inPlace(Sample* arr) {
            Sample sum = 0;
            for (int i = 0; i < size; ++i) {
                sum += arr[i];
            }
            sum *= multiplier;
            for (int i = 0; i < size; ++i) {
                arr[i] += sum;
            }
        }
    };

    // Use like `Hadamard<double, 8>::inPlace(data)` - size must be a power of 2
    template <typename Sample, int size>
    class [[maybe_unused]] Hadamard {
    public:
        static inline void recursiveUnscaled(Sample* data) {
            if constexpr (size <= 1)
                return;
            else {
                constexpr int hSize = size / 2;
                // Two (unscaled) Hadamards of half the size
                Hadamard<Sample, hSize>::recursiveUnscaled(data);
                Hadamard<Sample, hSize>::recursiveUnscaled(data + hSize);
                // Combine the two halves using sum/difference
                for (int i = 0; i < hSize; ++i) {
                    Sample a = data[i];
                    Sample b = data[i + hSize];
                    data[i] = (a + b);
                    data[i + hSize] = (a - b);
                }
            }
        }
        static inline void inPlace(Sample* data) {
            recursiveUnscaled(data);
            auto scalingFactor = static_cast<Sample>(std::sqrt(1.0 / size));
            for (int c = 0; c < size; ++c) {
                data[c] *= scalingFactor;
            }
        }
    };
} // namespace SDSP