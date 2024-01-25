//
// Created by Syl Morrison on 14/04/2023.
//
#include <catch2/catch_test_macros.hpp>
#include "../Types/CircularBuffer.h"
#include <juce_core/juce_core.h>

namespace SDSP::Testing {
    bool testCircularBufferIR(int delay) {
        SDSP::CircularBuffer<float> buffer{ 0.0, 44100 };
        buffer.prepare(256, 44100.0);
        buffer.setDelay(delay);
        std::vector<float> impulse{ 1.0f };
        for (; delay > 0; --delay) {
            impulse.emplace_back(0.0f);
        }
        std::vector<float> res;
        res.reserve(impulse.size());
        for (auto& s : impulse) {
            res.emplace_back(buffer.getNextSample(s));
        }
        return res.back() == 1;
    }

    TEST_CASE("Test Circular Buffer", "[circularBuffer]") {
        SECTION("Test out of range delay times") {
            CircularBuffer<float> buffer(0, 0.0);
            buffer.prepare(256, 44100.0);
            buffer.setMaxDelaySeconds(1);
            buffer.setDelay(std::numeric_limits<int>::max());
            REQUIRE(buffer.getDelay() == 44100);
            REQUIRE(buffer.getNextSample(0) == 0);
            buffer.setDelay(std::numeric_limits<int>::min());
            REQUIRE(buffer.getDelay() == 0);
            REQUIRE(buffer.getNextSample(0) == 0);
        }
        // TODO: Figure out how to test interpolating delay times

        SECTION("Sanity check buffer output") {
            REQUIRE(testCircularBufferIR(10));
            REQUIRE(testCircularBufferIR(100));
            REQUIRE(testCircularBufferIR(512));
        }
    }
} // namespace SDSP::Testing
