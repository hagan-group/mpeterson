#pragma once

template <class Generator>
class RandomForce {
   public:
    RandomForce(Generator rng) : rng_{rng} {}

    Scalar energy() const { return 0.0; }

    Vector force() { return {rng_(), rng_(), rng_()}; }

   private:
    Generator rng_;
};