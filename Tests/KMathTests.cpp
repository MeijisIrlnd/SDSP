//
// Created by Syl Morrison on 14/04/2023.
//
#include "../KMath.h"
#include <catch2/catch_test_macros.hpp>
namespace SDSP::Testing {
    TEST_CASE( "Test SDSP::KMath::Lerp<float>", "[lerp]") {
        REQUIRE(SDSP::KMath::Lerp<float>(0, 1, 0) == 0);
        REQUIRE(SDSP::KMath::Lerp<float>(0, 1, 1) == 1);
        REQUIRE(SDSP::KMath::Lerp<float>(0, 1, -1) == 0);
        REQUIRE(SDSP::KMath::Lerp<float>(0, 1, 2) == 1);
    }

    TEST_CASE( "Test SDSP::KMath::isPrime", "[isPrime]") {
        REQUIRE(SDSP::KMath::isPrime(-2) == false);
        REQUIRE(SDSP::KMath::isPrime(0) == false);
        REQUIRE(SDSP::KMath::isPrime(1) == false);
        REQUIRE(SDSP::KMath::isPrime(4) == false);
        REQUIRE(SDSP::KMath::isPrime(6) == false);
        for(auto& p : KMath::s_primes) {
            REQUIRE(SDSP::KMath::isPrime(p) == true);
        }
    }

    TEST_CASE("Test SDSP::KMath::nearestPrime", "[nearestPrime]") {
        REQUIRE(SDSP::KMath::nearestPrime(0) == 2);
        REQUIRE(SDSP::KMath::nearestPrime(1) == 2);
        REQUIRE(SDSP::KMath::nearestPrime(-1) == 2);
        REQUIRE(SDSP::KMath::nearestPrime(4) == 3);
        REQUIRE(SDSP::KMath::nearestPrime(530) == 523);
    }

    TEST_CASE("Test SDSP::KMath::Log_n", "[log]") {
        REQUIRE_THROWS(SDSP::KMath::log(0, 2));
        REQUIRE(SDSP::KMath::log(1, 0) == 0);
        REQUIRE(SDSP::KMath::log(2.0f, 0.5f) == -1);
        REQUIRE(SDSP::KMath::log(64, 2) == std::log2(64));
    }

    TEST_CASE("Test SDSP::getNearestCoprime", "[getNearestCoprime]") {
        REQUIRE(SDSP::KMath::getNearestCoprime(0) == 2);
        REQUIRE(SDSP::KMath::getNearestCoprime(54) == 53);

    }
}