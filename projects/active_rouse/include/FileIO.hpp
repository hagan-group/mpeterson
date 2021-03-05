/**
 * Utilities for reading and writing CSV files
 */

#pragma once

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <type_traits>

// Generic type cpnvertion is simply a static_cast
template <class To, class From>
To to(const From& from) {
    if constexpr (std::is_same_v<To, std::string>) {
        return std::to_string(from);
    } else if constexpr (std::is_same_v<From, std::string>) {
        auto out = To{};
        std::stringstream stream{from};
        stream >> out;
        return out;
    } else {
        return static_cast<To>(from);
    }
}

// CSV parsing, optionally with conversion of cells to a particular type
template <class CellType = std::string>
auto read_csv(const std::string& infile) {
    using Row = std::vector<CellType>;
    using CsvOutput = std::vector<Row>;

    std::ifstream in(infile);
    if (!in.good()) {
        throw std::runtime_error("Unable to read file " + infile);
    }

    CsvOutput output;

    std::string line;
    std::string cell;
    std::stringstream line_stream;

    size_t row = 0;
    while (std::getline(in, line)) {
        output.push_back(Row{});

        line_stream = std::stringstream(line);
        while (std::getline(line_stream, cell, ',')) {
            output[row].push_back(to<CellType>(cell));
        }

        ++row;
    }

    return output;
}

// Writing of CSV from any doubly nested container
template <class Container>
void write_csv(std::string fname, const Container& cont) {
    std::ofstream out(fname);

    for (const auto& row : cont) {
        for (const auto& cell : row) {
            out << cell << ",";
        }
        out << "\n";
    }
}
