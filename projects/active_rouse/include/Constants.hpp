#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"

#include <limits>

// scalars
constexpr Scalar Pi = 3.141592653589793;
constexpr Scalar E = 2.718281828459045;
constexpr Scalar Inf = std::numeric_limits<double>::infinity();
constexpr Scalar NaN = std::numeric_limits<double>::quiet_NaN();

// identity matrix
constexpr Matrix Identity = {1, 0, 0, 
                             0, 1, 0,
                             0, 0, 1};
