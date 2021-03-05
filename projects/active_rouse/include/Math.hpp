#pragma once

#include "Constants.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

#include <cmath>
#include <tuple>

// ==============================================
// Dot products (v-v, v-m, m-v, m-m)
// ==============================================

constexpr Vector::Element dot(const Vector& u, const Vector& v) noexcept {
    return u(0) * v(0) + u(1) * v(1) + u(2) * v(2);
}

Vector dot(const Matrix& m, const Vector& v) noexcept {
    Vector out;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            out(i) += m(i, j) * v(j);
        }
    }
    return out;
}

Vector dot(const Vector& v, const Matrix& m) noexcept {
    Vector out;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            out(j) += v(i) * m(i, j);
        }
    }
    return out;
}

Matrix dot(const Matrix& a, const Matrix& b) noexcept {
    Matrix c;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return c;
}

// ==============================================
// Out products (v-v only)
// ==============================================

Matrix outer(const Vector& u, const Vector& v) noexcept {
    Matrix out;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            out(i, j) = u(i) * v(j);
        }
    }
    return out;
}

// ==============================================
// Vector functionality
// ==============================================

Vector normalize(Vector u) {
    u.normalize();
    return u;
}

Scalar angle(const Vector& u, const Vector& v) noexcept {
    auto d = dot(u, v);
    if (d == 0.0) return 0.5 * Pi;

    d /= std::sqrt(u.norm2() * v.norm2());
    if (d >= 1.0) return 0.0;
    if (d <=-1.0) return Pi;
    return std::acos(d);
}

// ==============================================
// Matrix functionality
// ==============================================

constexpr Matrix::Element trace(const Matrix& m) noexcept {
    return m(0, 0) + m(1, 1) + m(2, 2);
}

constexpr Matrix::Element det(const Matrix& m) noexcept {
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
           m(1, 1) * (m(0, 0) * m(2, 2) - m(2, 0) * m(0, 2)) +
           m(2, 2) * (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1));
}

// Construct matrix from row vectors
constexpr Matrix as_rows(const Vector& u, const Vector& v,
                         const Vector& w) noexcept {
    return {u(0), u(1), u(2), v(0), v(1), v(2), w(0), w(1), w(2)};
    // return out;
}

// Construct matrix from column vectors
constexpr Matrix as_cols(const Vector& u, const Vector& v,
                         const Vector& w) noexcept {
    return {u(0), v(0), w(0), u(1), v(1), w(1), u(2), v(2), w(2)};
    // return out;
}

constexpr auto from_rows(const Matrix& m) noexcept {
    Vector u = {m(0, 0), m(0, 1), m(0, 2)};
    Vector v = {m(1, 0), m(1, 1), m(1, 2)};
    Vector w = {m(2, 0), m(2, 1), m(2, 2)};
    return std::make_tuple(u, v, w);
}

constexpr auto from_cols(const Matrix& m) noexcept {
    Vector u = {m(0, 0), m(1, 0), m(2, 0)};
    Vector v = {m(0, 1), m(1, 1), m(2, 1)};
    Vector w = {m(0, 2), m(1, 2), m(2, 2)};
    return std::make_tuple(u, v, w);
}

// ==============================================
// Miscellaneous
// ==============================================

template <class T>
constexpr int kronecker_delta(const T& x, const T& y) noexcept {
    return x == y ? 1 : 0;
}

template <class T>
constexpr int sign(const T& x) noexcept {
    return (x > T()) - (x < T());
}
