
//
// Created by Syl on 17/04/2023.
//
#include <catch2/catch_test_macros.hpp>
#include "../Helpers/Helpers.h"
namespace SDSP::Testing
{
    TEST_CASE("Test zero_array", "[zero_array]") {
        std::array<float, 2> arr1{ 1, 1};
        Helpers::zero_array(arr1);
        REQUIRE((arr1[0] == 0 && arr1[1] == 0));
    }

    TEST_CASE("Test copy_array", "[copy_array]") {
        std::array<float, 2> arr1{1, 1}, arr2{0, 0};
        Helpers::copy_array(arr1, arr2);
        REQUIRE((arr1[0] == 0 && arr1[1] == 0));
    }

    TEST_CASE("Test fill_array", "[fill_array]") {
        std::array<float, 2> arr{0, 0};
        Helpers::fill_array(arr, 1);
        REQUIRE((arr[0] == 1 && arr[1] == 1));
    }



}