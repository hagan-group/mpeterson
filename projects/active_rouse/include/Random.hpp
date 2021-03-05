#pragma once

#include <random>

template <class Distribution>
class Random {
   public:
    Random() : seed_{std::random_device()()}, gen_{seed_} {}

    Random(std::uint64_t seed) : seed_{seed}, gen_{seed_} {}

    Random(const Random&) = default;
    Random(Random&&) = default;

    auto seed() const { return seed_; }

    double operator()() { return dist_(gen_); }

   private:
    const std::uint64_t seed_;
    std::mt19937_64 gen_;
    Distribution dist_;
};


using RandomNormal = Random<std::normal_distribution<double>>;
using RandomUniform = Random<std::uniform_real_distribution<double>>;
