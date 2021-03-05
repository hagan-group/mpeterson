#pragma once

#include "Math.hpp"
// #include "NeighborList.hpp"

#include <tuple>
#include <vector>

struct Atoms {
    using PairList = List<std::tuple<size_t, size_t>>;
    using BondList = List<std::tuple<size_t, size_t>>;
    using AngleList = List<std::tuple<size_t, size_t, size_t>>;
    
    // data members
    List<Vector> pos;
    List<Vector> vel;
    List<Vector> force;
    List<size_t> mol;
    PairList pairs;
    BondList bonds;
    AngleList angles;

    // apparently we need to explicitly default these?
    Atoms() = default;
    Atoms(const Atoms&) = default;
    Atoms(Atoms&&) = default;
    Atoms& operator=(const Atoms&) = default;
    Atoms& operator=(Atoms&&) = default;

    // size constructor
    Atoms(size_t size) : pos(size), vel(size), force(size), mol(size) {}

    size_t size() const { return pos.size(); }

    void add_atom(Vector r = Vector{}) {
        pos.push_back(r);
        // vel.push_back(v);
        force.push_back(Vector{});
    }
};
