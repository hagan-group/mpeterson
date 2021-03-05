#pragma once

#include "Constants.hpp"
#include "Math.hpp"
#include "Types.hpp"

#include <cmath>

class HarmonicBond {
   public:
    using Result = std::tuple<Vector, Vector>;

    HarmonicBond(Scalar k = 1.0, Scalar b = 0.0) : k_{k}, b_{b} {}

    Scalar energy(const Vector& u, const Vector& v) const noexcept {
        return 0.5 * k_ * std::pow(((v - u).norm() - b_), 2);
    }

    Result force(const Vector& u, const Vector& v) const noexcept {
        auto rhat = v - u;
        auto r = rhat.normalize();
        auto out = k_ * (r - b_) * rhat;
        return {out, -out};
    }

   private:
    Scalar k_;
    Scalar b_;
};
