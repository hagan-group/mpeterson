#pragma once

#include "Constants.hpp"
#include "Math.hpp"
#include "Types.hpp"

#include <cmath>
#include <tuple>

// a Functor that computes the Lennard-Jones potential (or force) as a function
// of the position of two atoms
class LennardJones {
   public:
    using Result = std::tuple<Vector, Vector>;

    LennardJones(Scalar epsilon = 1.0, Scalar sigma = 1.0, Scalar cutoff = 2.5,
                 bool shift = true)
        : a_{4.0 * epsilon * std::pow(sigma, 12)},
          b_{4.0 * epsilon * std::pow(sigma, 6)},
          cutoff2_{cutoff * cutoff},
          shift_{0.0} {
        if (shift) {
            auto r6 = std::pow(cutoff, 6);
            auto r12 = r6 * r6;
            shift_ = a_ / r12 - b_ / r6;
        }
    }

    Scalar energy(const Vector& u, const Vector& v) const noexcept {
        auto r2 = (v - u).norm2();
        if (r2 > 0.0 && r2 < cutoff2_) {
            auto r6 = std::pow((v - u).norm2(), 3);
            auto r12 = r6 * r6;
            return a_ / r12 - b_ / r6 - shift_;
        } else if (r2 == 0.0) {
            return Inf;
        } else {
            return 0.0;
        }
    }

    Result force(const Vector& u, const Vector& v) const noexcept {
        auto rvec = v - u;
        auto r2 = rvec.norm2();
        if (r2 > 0.0 && r2 < cutoff2_) {
            auto r8 = std::pow(r2, 4);
            auto r14 = r8 * std::pow(r2, 3);
            auto out = rvec * (6 * b_ / r8 - 12 * a_ / r14);
            return {out, -out};
        } else {
            return {};
        }
    }

   private:
    Scalar a_;
    Scalar b_;
    Scalar cutoff2_;
    Scalar shift_;
};
