//
// Created by Syl on 18/04/2023.
//
#include "../Oscillators/SDSPOscillator.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Test SDSPOscillator", "[SDSPOscillator]") {
    double sampleRate{ 100 };
    SECTION("Test Generators") {
        SECTION("Test Sine") {
            SDSP::Oscillators::SDSPOscillator oscillator{ false, SDSP::Oscillators::SHAPE::SINE };
            oscillator.prepareToPlay(256, sampleRate);
            oscillator.setFrequency(1.0f);
            auto phaseIncrement = 1 / static_cast<float>(sampleRate);
            auto phase = 0.0f;
            for (auto i = 0; i < 256; ++i) {
                REQUIRE_THAT(oscillator.processSample(),
                             Catch::Matchers::WithinRel(std::sin(phase * juce::MathConstants<float>::twoPi)));
                phase += phaseIncrement;
                phase = std::fmod(phase, 1.0f);
                REQUIRE_THAT(oscillator.getPhase(), Catch::Matchers::WithinRel(phase));
            }
        }
        SECTION("Test Square") {
            SDSP::Oscillators::SDSPOscillator oscillator(false, SDSP::Oscillators::SHAPE::SQUARE);
            oscillator.prepareToPlay(256, sampleRate);
            oscillator.setFrequency(1.0f);
            // Well okay, with a frequency of 1hz, we should have sr / 2 1s, and sr/2 -1s..
            std::vector<float> testSamples(static_cast<size_t>(sampleRate));
            std::fill(testSamples.begin(), testSamples.begin() + static_cast<long long>(sampleRate) / 2 + 1, -1.0f);
            std::fill(testSamples.begin() + static_cast<long long>(sampleRate) / 2 + 1, testSamples.end(), 1.0f);
            std::vector<float> res;
            for (float testSample : testSamples) {
                REQUIRE_THAT(oscillator.processSample(), Catch::Matchers::WithinRel(testSample));
            }
        }
        SECTION("Test Triangle") {
            SDSP::Oscillators::SDSPOscillator oscillator(false, SDSP::Oscillators::SHAPE::TRI);
            oscillator.prepareToPlay(256, sampleRate);
            oscillator.setFrequency(1.0f);
            auto phaseIncrement = 1.0f / static_cast<float>(sampleRate);
            auto phase = 0.0f;
            auto triGen = [](float x) {
                auto numerator = std::asin(std::sin(x * juce::MathConstants<float>::twoPi));
                return numerator / juce::MathConstants<float>::halfPi;
            };
            for (auto sample = 0; sample < 256; ++sample) {
                REQUIRE_THAT(oscillator.processSample(), Catch::Matchers::WithinRel(triGen(phase)));
                phase += phaseIncrement;
                phase = std::fmod(phase, 1.0f);
            }
        }
        SECTION("Test Saw") {
            SDSP::Oscillators::SDSPOscillator oscillator(false, SDSP::Oscillators::SHAPE::SAW);
            oscillator.prepareToPlay(256, sampleRate);
            oscillator.setFrequency(1.0f);
            auto phaseIncrement = 1.0f / static_cast<float>(sampleRate);
            auto phase = 0.0f;
            for (auto sample = 0; sample < 256; ++sample) {
                REQUIRE_THAT(oscillator.processSample(), Catch::Matchers::WithinRel((2 * phase) - 1));
                phase += phaseIncrement;
                phase = std::fmod(phase, 1.0f);
            }
        }
    }
    SECTION("Test params") {
        SECTION("Test phase") {
            std::vector<float> testFreqs{ 1.0f, 1000.0f, 250.0f };
            SDSP::Oscillators::SDSPOscillator oscillator(false, SDSP::Oscillators::SHAPE::SINE);
            oscillator.prepareToPlay(256, sampleRate);
            for (auto& f : testFreqs) {
                oscillator.retrigger();
                REQUIRE(oscillator.getPhase() == 0.0f);
                oscillator.setFrequency(f);
                auto phaseIncrement = f / static_cast<float>(sampleRate);
                auto phase = 0.0f;
                for (auto i = 0; i < 256; ++i) {
                    REQUIRE_THAT(oscillator.getPhase(), Catch::Matchers::WithinRel(phase));
                    phase += phaseIncrement;
                    phase = std::fmod(phase, 1.0f);
                    [[maybe_unused]] auto _ = oscillator.processSample();
                }
            }
        }
    }
}