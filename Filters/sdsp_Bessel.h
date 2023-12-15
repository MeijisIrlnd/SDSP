#ifndef SDSP_BESSEL_H
#define SDSP_BESSEL_H
#include "DSPBiquad.h"
#include "SmoothedFilterCoefficients.h"
#include "RBJCoefficients.h"
#include <vector>
namespace SDSP::Filters {
    struct BesselFilterConfig {
        float frequencyMultiplier, q;
    };
    static const std::array<std::vector<BesselFilterConfig>, 4> s_besselFilterConfigs = {
        std::vector<BesselFilterConfig>{ BesselFilterConfig{ 1.27201964951f, 0.57735026919 } },
        std::vector<BesselFilterConfig>{ BesselFilterConfig{ 1.60335751622f, 0.805538281842f }, BesselFilterConfig{ 1.43017155999f, 0.521934581669f } },
        std::vector<BesselFilterConfig>{ BesselFilterConfig{ 1.9047076123f, 1.02331395383f }, BesselFilterConfig{ 1.68916826762f, 0.611194546878f }, BesselFilterConfig{ 1.60391912877f, 0.510317824749f } },
        std::vector<BesselFilterConfig>{ BesselFilterConfig{ 2.18872623053f, 1.22566942541f }, BesselFilterConfig{ 1.95319575902f, 0.710852074442f }, BesselFilterConfig{ 1.8320926012f, 0.559609164796f }, BesselFilterConfig{ 1.77846591177f, 0.505991069397f } }
    };

    template <size_t Order>
    static std::vector<BesselFilterConfig> getBesselConfigsForOrder() {
        static_assert(Order % 2 == 0 && Order <= 8);
        size_t actualIndex = Order / 2;
        return s_besselFilterConfigs[actualIndex];
    }

    enum class BESSEL_FILTER_TYPE {
        LOWPASS,
        LOW_SHELF,
        BELL,
        HIGH_SHELF,
        HIGHPASS
    };
    template <size_t Order, BESSEL_FILTER_TYPE FilterType>
    class Bessel {
    public:
        struct BesselParameters {
            float cf;
            float gain;
        };
        Bessel() : m_filterParams(getBesselConfigsForOrder<Order>()) {
            //
        }

        void prepareToPlay(int /* samplesPerBlockExpected */, double sampleRate, BesselParameters&& defaultParams) {
            m_sampleRate = sampleRate;
            checkFilters(std::move(defaultParams));
            for (auto stage = 0; stage < static_cast<int>(Order / 2); ++stage) {
                std::memcpy(m_coeffs.current(stage), m_coeffs.target(stage), sizeof(double) * 6);
                m_filter.setCoefficients(m_coeffs.target(stage), stage);
            }
            m_samplesUntilFilterUpdate = 0;
        }

        [[nodiscard]] float processSample(float x, BesselParameters&& params) {
            checkFilters(std::move(params));
            const auto filtered = m_filter.processSample(x);
            return filtered;
        }

    private:
        void checkFilters(BesselParameters&& params) {
            if (m_samplesUntilFilterUpdate == 0) {
                setFilterTargets(std::move(params));
                m_samplesUntilFilterUpdate = m_filterUpdateRate;
            }
            m_coeffs.interpolate();
            for (size_t stage = 0; stage < Order / 2; ++stage) {
                const auto iStage = static_cast<int>(stage);
                m_filter.setCoefficients(m_coeffs.current(iStage), iStage);
            }
            --m_samplesUntilFilterUpdate;
        }

        void setFilterTargets(BesselParameters&& params) noexcept {
            for (size_t i = 0; i < Order / 2; ++i) {
                const auto [multiplier, q] = m_filterParams[i];
                auto scaledCf = params.cf * multiplier;
                // problem here - if scaledCf goes above nyquist, we're fucked. We could either oversample, OR:
                scaledCf = scaledCf > 20000.0f ? 20000.0f : scaledCf;
                const auto iStage{ static_cast<int>(i) };
                if constexpr (FilterType == BESSEL_FILTER_TYPE::LOWPASS) {
                    SDSP::RBJ::lowpass(m_coeffs.target(iStage), m_sampleRate, scaledCf, q);
                } else if constexpr (FilterType == BESSEL_FILTER_TYPE::LOW_SHELF) {
                    SDSP::RBJ::lowShelf(m_coeffs.target(iStage), m_sampleRate, scaledCf, params.gain, q);
                } else if constexpr (FilterType == BESSEL_FILTER_TYPE::BELL) {
                    SDSP::RBJ::bellWithQ(m_coeffs.target(iStage), m_sampleRate, params.gain, scaledCf, q);
                } else if constexpr (FilterType == BESSEL_FILTER_TYPE::HIGH_SHELF) {
                    SDSP::RBJ::highShelf(m_coeffs.target(iStage), m_sampleRate, scaledCf, params.gain, q);
                } else if constexpr (FilterType == BESSEL_FILTER_TYPE::HIGHPASS) {
                    SDSP::RBJ::highpass(m_coeffs.target(iStage), m_sampleRate, scaledCf, q);
                }
            }
        }

        double m_sampleRate{};
        const int m_filterUpdateRate{ 100 };
        int m_samplesUntilFilterUpdate{ 0 };
        const std::vector<BesselFilterConfig> m_filterParams;
        SDSP::BiquadCascade<Order / 2> m_filter;
        SDSP::SmoothedFilterCoefficients<Order / 2> m_coeffs;
    };
} // namespace SDSP::Filters
#endif