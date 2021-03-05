#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

class DumpWriter {
   public:
    template <class AtomList>
    void write_csv(std::string filename, const AtomList& atoms) {
        auto file = std::ofstream(filename);
        
        file << std::fixed;
        file << std::setprecision(5);
        for (const auto& r : atoms.pos) {
            file << r(0) << "," << r(1) << "," << r(2) << "\n";
        }
    }
};
