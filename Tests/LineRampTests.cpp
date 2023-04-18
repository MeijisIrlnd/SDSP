//
// Created by Syl on 17/04/2023.
//
#include "../Utils/LineRamp.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
namespace SDSP::Testing {
    TEST_CASE("Test LineRamp", "[LineRamp]") {
        SECTION("Test t = 0") {
            LineRamp<float> ramp;
            ramp.prepare(44100);
            ramp.set(1,0);
            for(auto i = 0; i < 10; ++i) {
                REQUIRE(ramp.process() == 1);
            }
        }
        SECTION("Test invalid times") {
            LineRamp<float> ramp;
            ramp.prepare(44100);
            ramp.set(1, -1);
            for(auto i = 0; i < 10; ++i) {
                REQUIRE(ramp.process() == 1);
            }
        }

        SECTION("Test interpolation values") {
            LineRamp<float> ramp;
            ramp.prepare(44100);
            ramp.set(0, 1, 1);
            auto nSamplesToProcess = 44100 / 1000;
            auto increment = 1 / static_cast<float>(nSamplesToProcess);
            for (auto i = 0; i < nSamplesToProcess; ++i) {
                REQUIRE_THAT(ramp.process(), Catch::Matchers::WithinRel(increment * i));
            }
        }

        SECTION("Test stop") {
            LineRamp<float> ramp;
            ramp.prepare(44100);
            ramp.set(0, 1, 1);
            auto nSamplesToProcess = 44100 / 1000;
            auto increment = 1 / static_cast<float>(nSamplesToProcess);
            for(auto i = 0; i < nSamplesToProcess / 2; ++i) {
                ramp.process();
            }
            ramp.stop();
            for(auto i = 0; i < 10; ++i) {
                REQUIRE_THAT(ramp.process(), Catch::Matchers::WithinRel(increment * (nSamplesToProcess / 2)));
            }
        }
    }
}