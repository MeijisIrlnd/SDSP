//
// Created by Syl Morrison on 28/07/2023.
//

#ifndef COLLAGE_SDSPWAVETABLE_H
#define COLLAGE_SDSPWAVETABLE_H
#include "OscillatorShape.h"
#include "SDSPOscillator.h"
#include <functional>
namespace SDSP::Oscillators {
    template<size_t NumPoints, SHAPE Shape>
    class [[maybe_unused]] SDSPWavetable {
    public:

        template<typename T = void, typename = std::enable_if_t<Shape != SHAPE::CUSTOM, T> >
        SDSPWavetable() {
            generateWavetable();
        }

        template<typename T = void, typename = std::enable_if_t<Shape == SHAPE::CUSTOM, T> >
        explicit SDSPWavetable(const std::function<float(float)>& generatorFunction) : m_generatorFunction(generatorFunction) {
            generateWavetable();
        }

        void prepareToPlay(int /*samplesPerBlockExpected*/, double sampleRate) {
            m_sampleRate = sampleRate;
        }

        [[nodiscard]] float processSample() noexcept {
            // Aw shit here we go again
            auto idx0 = static_cast<size_t>(m_currentIndex);
            auto idx1 = idx0 + 1;
            auto frac = m_currentIndex - static_cast<float>(idx0);
            auto v0 = m_table[idx0];
            auto v1 = m_table[idx1];
            auto interpolated = v0 + frac * (v1 - v0);
            if((m_currentIndex += m_tableDelta) > static_cast<float>(NumPoints)) {
                m_currentIndex -= static_cast<float>(NumPoints);
            }
            return interpolated;
        }

        SDSP_INLINE void setFrequency(float newFrequency) {
            auto tsOverSr = static_cast<float>(NumPoints) / static_cast<float>(m_sampleRate);
            m_tableDelta = newFrequency * tsOverSr;
        }


    private:
        void generateWavetable() {
            SDSPOscillator osc{ false };
            osc.setShape(Shape);
            osc.prepareToPlay(512, static_cast<double>(NumPoints));
            osc.setFrequency(1.0f);
            if constexpr(Shape == SHAPE::CUSTOM) {
                osc.setFunction(m_generatorFunction);
            }
            for(size_t i = 0; i < NumPoints; ++i) {
                m_table[i] = osc.processSample();
            }
            // Wraparound point
            m_table[NumPoints] = m_table[0];
        }

        double m_sampleRate{};
        std::function<float(float)> m_generatorFunction{ nullptr };
        std::array<float, NumPoints + 1> m_table;
        float m_currentIndex{ 0.0f}, m_tableDelta{ 0.0f };

    };
}
#endif //COLLAGE_SDSPWAVETABLE_H
