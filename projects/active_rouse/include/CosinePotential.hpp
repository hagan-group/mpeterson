#pragma once

#include "Math.hpp"
#include "Types.hpp"

#include <cmath>
#include <tuple>

/**
 * CosinePotential
 *
 * Given vectors u, v, and w, computes the energy/force based on a potential of
 * the form U = k * (1 - x.y) where x = normalized(v - u) and y = normalized(w -
 * v).
 */
class CosinePotential {
   public:
    using Result = std::tuple<Vector, Vector, Vector>;

    CosinePotential(Scalar k) : k_{k} {}

    Scalar energy(const Vector& u, const Vector& v, const Vector& w) const
        noexcept {
        return k_ * (1 - dot(normalize(v - u), normalize(w - v)));
    }

    Result force(const Vector& u, const Vector& v, const Vector& w) const
        noexcept {
        // vectors
        auto rij = v - u;
        auto rjk = w - v;

        // magnitudes
        auto mij = rij.normalize();
        auto mjk = rjk.normalize();

        auto Qij = Identity - outer(rij, rij);
        auto Qjk = Identity - outer(rjk, rjk);

        auto fi = mij == 0.0 ? Vector{} : -k_ * dot(Qij, rjk / mij);
        auto fk = mjk == 0.0 ? Vector{} :  k_ * dot(Qjk, rij / mjk);
        
        return {fi, -(fi + fk), fk};
    }

   private:
    Scalar k_;
};
