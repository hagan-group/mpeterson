#include <cassert>
#include <experimental/filesystem>
#include <iostream>
#include <random>
#include <algorithm>

#include "Atoms.hpp"
#include "Constants.hpp"
#include "DumpWriter.hpp"
#include "FileIO.hpp"
#include "Math.hpp"
#include "NeighborList.hpp"
#include "Potentials.hpp"
#include "Random.hpp"
#include "Types.hpp"
#include "Vector.hpp"

namespace fs = std::experimental::filesystem;

// ==========================================
// Simulation parameters
// ==========================================

// Enable features
constexpr bool ENABLE_EXPERIMENTAL = true;
constexpr bool CONVERT_ACTIVITY = true;

// System parameters
constexpr size_t NDIMS = 3;

// Atom parameters
constexpr size_t NATOMS = 51;
constexpr size_t NPOLYS = 50;

// RNG seed
constexpr size_t SEED = 0x397446F70D57E555UL;

// WCA parameters
constexpr double WCA_EPSILON = 1.0;
constexpr double WCA_SIGMA = 1.0;
constexpr double WCA_CUTOFF = 1.122462048309373;

// Bond parameters
constexpr double BOND_STRENGTH = 200;
constexpr double BOND_LENGTH = 1.0;

// Run parameters
// constexpr double TIMESTEP = 0.001;
// constexpr size_t NSTEPS = 100'000'000;
// constexpr size_t DUMPSTEP = 100'000;
constexpr double TIMESTEP = 0.0001;
constexpr size_t NSTEPS = 1'000'000'000;
constexpr size_t DUMPSTEP = 1'000'000;
constexpr bool CONTINUE = true;

// derived parameters
constexpr double NOISE_COEFF = std::sqrt(2.0 * TIMESTEP);
constexpr double NOISE_VAR = NDIMS * 2.0 * TIMESTEP;
constexpr double SYS_SKIN = WCA_CUTOFF;
constexpr double SYS_CUTOFF = WCA_CUTOFF + SYS_SKIN;

// Relevant directories
auto CWD = fs::path("./");

// Function declarations
void initialize(Atoms&, size_t&);
void build_neighbor_list(Atoms&);

// ==========================================
// Main
// ==========================================

int main(int argc, char* argv[]) {
    // ==========================================
    // Initialization
    // ==========================================

    // Get the active force strength and bend force coefficient from command
    // line arguments
    double activity = 0.0;
    double stiffness = 0.0;

    if (argc > 1) activity = std::atof(argv[1]);
    if (argc > 2) stiffness = std::atof(argv[2]);
    if (argc > 3) {
        std::cerr << "Usage: rouse ACTIVITY STIFFNESS\n";
        return 1;
    }

    Atoms atoms(NATOMS * NPOLYS);
    size_t step = 0;
    double dist2_moved = 0.0;
    double max_dist2 = 0.0;

    initialize(atoms, step);

    if constexpr (CONVERT_ACTIVITY) {
        // alpha = fa * N / 2 * k = fa * (N / 6)
        // --> fa = 6 * activity / N
        // should N be NATOMS or NBONDS = NATOMS - 1?
        activity *= 6.0 / (NATOMS - 1);
    }

    auto wca = LennardJones(WCA_EPSILON, WCA_SIGMA, WCA_CUTOFF);
    auto bond = HarmonicBond(BOND_STRENGTH, BOND_LENGTH);
    auto bend = CosinePotential(stiffness);
    auto active = ActiveForce(activity);
    auto thermal = RandomForce(RandomNormal(SEED));

    // ==========================================
    // Run it!
    // ==========================================

    std::chrono::high_resolution_clock clock;
    auto start = clock.now();

    [[maybe_unused]] auto writer = DumpWriter();

    // size_t step = 0;
    while (step <= NSTEPS) {
        if (step % DUMPSTEP == 0) {
            #ifdef NDEBUG
            std::cout << "Writing dump file for step " << step << "\n";
            auto outfile = std::string("pos_") + std::to_string(step) + ".csv";
            writer.write_csv(outfile, atoms);
            
            #else
            std::cout << "At step " << step << std::endl;
            
            #endif
        }

        // update neighbor list if necessary
        if constexpr (ENABLE_EXPERIMENTAL) {
            if (dist2_moved >= 0.25 * SYS_SKIN * SYS_SKIN) {
                #ifndef NDEBUG
                std::cout << "Updating neighbor list at step " << step << std::endl;
                std::cout << "\tMoved distance of " << dist2_moved << std::endl;
                #endif
                
                build_neighbor_list(atoms);
                dist2_moved = 0.0;

                #ifndef NDEBUG
                for (auto pair : atoms.pairs) {
                    std::cout << std::get<0>(pair) << " <---> " << std::get<1>(pair) << std::endl;
                }
                #endif
            }
        }

        // compute pairwise forces
        for (auto [i, j] : atoms.pairs) {
            auto [fi, fj] = wca.force(atoms.pos[i], atoms.pos[j]);
            atoms.force[i] += fi;
            atoms.force[j] += fj;
        }

        // compute bond + active forces
        for (auto [i, j] : atoms.bonds) {
            auto [fi, fj] = bond.force(atoms.pos[i], atoms.pos[j]);
            auto [ai, aj] = active.force(atoms.pos[i], atoms.pos[j]);
            atoms.force[i] += fi + ai;
            atoms.force[j] += fj + aj;
        }

        // compute angle forces
        for (auto [i, j, k] : atoms.angles) {
            auto [fi, fj, fk] =
                bend.force(atoms.pos[i], atoms.pos[j], atoms.pos[k]);
            atoms.force[i] += fi;
            atoms.force[j] += fj;
            atoms.force[k] += fk;
        }

        // update atom positions
        max_dist2 = 0.0;
        for (size_t i = 0; i < atoms.size(); ++i) {
            auto move = TIMESTEP * atoms.force[i] + NOISE_COEFF * thermal.force();
            atoms.pos[i] += move;
            atoms.force[i] = {};

            if (move.norm2() > max_dist2) {
                max_dist2 = move.norm2();
            }
        }

        dist2_moved += max_dist2;
        step += 1;
    }

    // print performance
    auto stop = clock.now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    auto steps_per_second = NSTEPS * 1e6 / time.count();

    std::cout << "Performance: " << steps_per_second << " steps/second\n";

    return 0;
}

// ==========================================
// Function definitions
// ==========================================

void initialize(Atoms& atoms, size_t& step) {
    if constexpr (CONTINUE) {
        std::cout << "Loading from previous run...\n";

        // first, determine the largest step currently available
        size_t max_step = 0;
        for (const auto& file : fs::directory_iterator(CWD)) {
            if (!fs::is_regular_file(file.path())) continue;
            if (!file.path().has_extension()) continue;

            auto name = file.path().filename().string();
            if (name.compare(0, 4, "pos_") != 0) continue;

            step = to<size_t>(name.substr(4, name.length() - 4));
            max_step = step > max_step ? step : max_step;
        }

        // Now load the corresponding file
        auto load_file =
            CWD / (std::string("pos_") + to<std::string>(max_step) +
                   std::string(".csv"));

        std::cout << "\tLoading file " << load_file << "\n";

        auto csv_data = read_csv<double>(load_file.string());
        for (size_t i = 0; i < NATOMS * NPOLYS; ++i) {
            atoms.pos[i].x() = csv_data[i][0];
            atoms.pos[i].y() = csv_data[i][1];
            atoms.pos[i].z() = csv_data[i][2];
        }

        step = max_step;
        csv_data.clear();
    } else {
        std::cout << "Initializing new run...\n";
        for (size_t p = 0; p < NPOLYS; ++p) {
            for (size_t i = 0; i < NATOMS; ++i) {
                size_t idx = NATOMS * p + i;

                atoms.pos[idx] = {to<double>(i), 0.0, 0.0};
            }
        }
    }

    for (size_t p = 0; p < NPOLYS; ++p) {
        for (size_t i = 0; i < NATOMS; ++i) {
            size_t idx = NATOMS * p + i;
            
            atoms.mol[idx] = p;
            
            for (size_t j = i + 2; j < NATOMS; ++j) {
                size_t idx2 = NATOMS * p + j;
                atoms.pairs.push_back({idx, idx2});
            }
            
            if (i < NATOMS - 1) atoms.bonds.push_back({idx, idx + 1});
            if (i < NATOMS - 2) atoms.angles.push_back({idx, idx + 1, idx + 2});
        }
    }

    // sort these, so that it's easier to search for elements later
    std::sort(atoms.bonds.begin(), atoms.bonds.end());
    std::sort(atoms.angles.begin(), atoms.angles.end());

    if constexpr (ENABLE_EXPERIMENTAL) {
        build_neighbor_list(atoms);
    }
}


void build_neighbor_list(Atoms& atoms) {

    Atoms::PairList pairs;

    Vector rij;
    std::tuple<size_t, size_t> pair;
    
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < i + NATOMS; ++j) {
            if (atoms.mol[i] != atoms.mol[j]) continue;

            rij = atoms.pos[j] - atoms.pos[i];
            if (rij.norm2() <= SYS_CUTOFF * SYS_CUTOFF) {
                pair = std::make_tuple(i, j);
                if (std::binary_search(atoms.bonds.begin(), atoms.bonds.end(), pair)) {
                    continue;
                }
                pairs.push_back(pair);
            }
        }
    }

    // swap atom pairlist
    atoms.pairs = std::move(pairs);
    
}
