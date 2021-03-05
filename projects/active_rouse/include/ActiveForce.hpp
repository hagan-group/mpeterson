#pragma once

#include "Constants.hpp"
#include "Math.hpp"
#include "Types.hpp"

#include <cmath>
#include <tuple>

class ActiveForce {
   public:
    using Result = std::tuple<Vector, Vector>;

    ActiveForce(Scalar alpha = 1.0) : alpha_{alpha} {}

    Scalar energy([[maybe_unused]] const Vector& u,
                  [[maybe_unused]] const Vector& v) {
        return 0.0;
    }

    Result force(const Vector& u, const Vector& v) noexcept {
        auto t = 0.5 * alpha_ * (v - u);
        // each atom attached to the bond feels the same force!
        // this is why an active force cannot be written as a simple
        // interaction potential --- the sum of forces is not 0!
        return {t, t};
    }

   private:
    Scalar alpha_;
};