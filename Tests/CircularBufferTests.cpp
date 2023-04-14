//
// Created by Syl Morrison on 14/04/2023.
//
#include <catch2/catch_test_macros.hpp>
#include "../Types/CircularBuffer.h"
namespace SDSP::Testing
{
    TEST_CASE("Test Circular Buffer", "[circularBuffer]") {
        SDSP::CircularBuffer<float> buffer{0, -1};

    }
}