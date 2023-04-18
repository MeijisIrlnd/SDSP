//
// Created by Syl on 18/04/2023.
//
#include "../Oscillators/SDSPOscillator.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Test SDSPOscillator", "[SDSPOscillator]") {
    SDSP::Oscillators::SDSPOscillator oscillator{ false, SDSP::Oscillators::SHAPE::SINE};
    oscillator.prepareToPlay(256, 44100);
    oscillator.setFrequency(1.0f);
    auto phaseIncrement = 1 / 44100.0f;
    auto phase = 0.0f;
    for(auto i = 0; i < 256; ++i) {
        REQUIRE_THAT(oscillator.processSample(), Catch::Matchers::WithinRel(std::sin(phase * juce::MathConstants<float>::twoPi)));
        phase += phaseIncrement;
        phase = std::fmod(phase, 1.0f);
        REQUIRE_THAT(oscillator.getPhase(), Catch::Matchers::WithinRel(phase));
    }

}