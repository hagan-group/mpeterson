#pragma once

#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>

// Alternate names for STL containers
template <class T, size_t N>
using Array = std::array<T, N>;

template <class T>
using List = std::vector<T>;

template <class T>
using Set = std::unordered_set<T>;

template <class Key, class Value>
using Dict = std::unordered_map<Key, Value>;

// basic mathematical types
using Scalar = double;
class Vector;
class Matrix;

// simulation types
struct Atoms;
class Domain;
class NeighborList;
