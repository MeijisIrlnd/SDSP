//
// Created by Syl on 17/04/2023.
//
#include "../Utils/LineRamp.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
namespace SDSP::Testing {
    TEST_CASE("Test LineRamp", "[LineRamp]") {
        LineRamp<float> ramp;
        ramp.prepare(44100);
        ramp.set(0,1, 1);
        auto nSamplesToProcess = 44100 / 1000;
        auto increment = 1 / static_cast<float>(nSamplesToProcess);
        for(auto i = 0; i < nSamplesToProcess; ++i) {
            REQUIRE_THAT(ramp.process(), Catch::Matchers::WithinRel(increment * i));
        }
    }
}