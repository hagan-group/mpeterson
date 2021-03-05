#include "Atoms.hpp"
#include "Constants.hpp"
#include "DumpWriter.hpp"
#include "Math.hpp"
#include "NeighborList.hpp"
#include "Potentials.hpp"
#include "Random.hpp"
#include "Vector.hpp"

#include <cassert>
#include <iostream>

int main() {
    // basic test setup will consist of 3 atoms making an "L" shape
    Atoms atoms(3);
    atoms.pos[0] = {0, 2, 0};
    atoms.pos[1] = {0, 0, 0};
    atoms.pos[2] = {1, 0, 0};

    // 3 pairwise interactions
    atoms.pairs.push_back({0, 1});
    atoms.pairs.push_back({0, 2});
    atoms.pairs.push_back({1, 2});

    // 2 bonds
    atoms.bonds.push_back({0, 1});
    atoms.bonds.push_back({1, 2});

    // 1 angle
    atoms.angles.push_back({0, 1, 2});

    // Basic forces
    auto lj = LennardJones(1.0, 1.0, 2.5);
    auto bend = CosinePotential(1.0);
    auto active = ActiveForce(1.0);
    auto bond = HarmonicBond(1.0, 1.5);
    auto random = RandomForce(RandomNormal(0x14354653));

    // ==========================================
    // Lennard-Jones Test
    // ==========================================

    auto [f0, f1] = lj.force(atoms.pos[0], atoms.pos[1]);
    Vector correct = {0.0, 0.181640625, 0.0};

    assert(f0 == -correct);
    assert(f1 == correct);

    std::tie(f0, f1) = lj.force(atoms.pos[1], atoms.pos[2]);
    correct = {24.0, 0.0, 0.0};

    assert(f0 == -correct);
    assert(f1 == correct);

    std::tie(f0, f1) = lj.force(atoms.pos[0], atoms.pos[2]);
    correct = {-0.0377856, 0.0755712, 0.0};

    assert(f0 == -correct);
    assert(f1 == correct);

    // ==========================================
    // Cosine Potential Test
    // ==========================================

    Vector f2;
    std::tie(f0, f1, f2) = bend.force(atoms.pos[0], atoms.pos[1], atoms.pos[2]);

    auto f0_c = Vector{-0.5, 0.0, 0.0};
    auto f2_c = Vector{0.0, -1.0, 0.0};
    auto f1_c = -(f0_c + f2_c);

    assert(f0 == f0_c);
    assert(f1 == f1_c);
    assert(f2 == f2_c);

    // ==========================================
    // Active Force Test
    // ==========================================

    std::tie(f0, f1) = active.force(atoms.pos[0], atoms.pos[1]);

    f0_c = {0.0, -1.0, 0.0};
    f1_c = f0_c;
    
    assert(f0 == f0_c);
    assert(f1 == f1_c);

    std::tie(f1, f2) = active.force(atoms.pos[1], atoms.pos[2]);

    f1_c = {0.5, 0.0, 0.0};
    f2_c = f1_c;

    assert(f1 == f1_c);
    assert(f2 == f2_c);

    
    // ==========================================
    // Bond Force Test
    // ==========================================

    std::tie(f0, f1) = bond.force(atoms.pos[0], atoms.pos[1]);

    f0_c = {0.0, -0.5, 0.0};
    f1_c = {0.0, +0.5, 0.0};

    assert(f0 == f0_c);
    assert(f1 == f1_c);

    std::tie(f1, f2) = bond.force(atoms.pos[1], atoms.pos[2]);

    f1_c = {-0.5, 0.0, 0.0};
    f2_c = {+0.5, 0.0, 0.0};

    assert(f1 == f1_c);
    assert(f2 == f2_c);

    std::cout << "All tests passed!\n";

    // ==========================================
    // Random Force Test
    // ==========================================

    // how to properly test this???
    List<Vector> vlist(2'000'000);

    for (int i = 0; i < 2'000'000; ++i) {
        vlist[i] = random.force();
    }

    Vector accum;
    for (auto v : vlist) {
        accum += v;
    }

    Vector mean = accum / 2'000'000;
    Vector std;
    accum = {0, 0, 0};
    for (auto v : vlist) {
        accum += (v - mean) * (v - mean);
    }
    std = accum / 2'000'000;

    std::cout << mean << std::endl;
    std::cout << std << std::endl;

}