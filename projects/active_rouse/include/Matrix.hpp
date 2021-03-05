#pragma once

#include "Types.hpp"

#include <array>
#include <cmath>
#include <cstdint>

using std::size_t;

#ifndef RELATIVE_TOL
#define RELATIVE_TOL 1e-10
#endif

#ifndef ABSOLUTE_TOL
#define ABSOLUTE_TOL 0.0
#endif

#define GENERATE_MATRIX_OP(OP)                                           \
    constexpr Matrix& operator OP##=(const Matrix& rhs) {                \
        for (size_t i = 0; i < c_.size(); ++i) {                         \
            c_[i] OP## = rhs.c_[i];                                      \
        }                                                                \
        return *this;                                                    \
    }                                                                    \
    constexpr Matrix& operator OP##=(Element rhs) {                      \
        for (size_t i = 0; i < c_.size(); ++i) {                         \
            c_[i] OP## = rhs;                                            \
        }                                                                \
        return *this;                                                    \
    }                                                                    \
    friend constexpr Matrix operator OP(Matrix lhs, const Matrix& rhs) { \
        lhs OP## = rhs;                                                  \
        return lhs;                                                      \
    }                                                                    \
    friend constexpr Matrix operator OP(Matrix lhs, Element rhs) {       \
        lhs OP## = rhs;                                                  \
        return lhs;                                                      \
    }                                                                    \
    friend constexpr Matrix operator OP(Element lhs, Matrix rhs) {       \
        for (size_t i = 0; i < rhs.c_.size(); ++i) {                     \
            rhs.c_[i] = lhs OP rhs.c_[i];                                \
        }                                                                \
        return rhs;                                                      \
    }

/// Matrix, class representing 3x3 matrices of real numbers
class Matrix {
   public:
    using Element = Scalar;
    using Container = Array<Element, 9>;

    // explicitly default normal constructors/assignment
    constexpr Matrix() noexcept : c_{} {}
    Matrix(const Matrix&) = default;
    Matrix(Matrix&&) = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&&) = default;

    template <class... Values,
              typename = std::enable_if_t<sizeof...(Values) == 9>>
    constexpr Matrix(Values... vals) noexcept : c_{Element(vals)...} {}

    // ==========================================
    // Utility
    // ==========================================

    constexpr size_t nrows() const noexcept { return 3; }
    constexpr size_t ncols() const noexcept { return 3; }
    constexpr size_t size() const noexcept { return 9; }

    // ==========================================
    // Mathematical functions
    // ==========================================

    constexpr Element norm2() const noexcept {
        Element out = 0.0;
        for (size_t i = 0; i < c_.size(); ++i) {
            out += c_[i] * c_[i];
        }
        return out;
    }

    Element norm() const noexcept { return std::sqrt(norm2()); }

    constexpr Matrix t() const noexcept {
        return {c_[0], c_[3], c_[6], c_[1], c_[4], c_[7], c_[2], c_[5], c_[8]};
    }

    // allow for elementwise mapping of a function
    template <class Function>
    constexpr Matrix& map(Function f) noexcept(noexcept(f(Element()))) {
        for (auto& v : c_) {
            v = f(v);
        }
        return *this;
    }
    template <class Function>
    constexpr Matrix map(Function f) const noexcept(noexcept(f(Element()))) {
        Matrix out;
        for (size_t i = 0; i < c_.size(); ++i) {
            out.c_[i] = f(c_[i]);
        }
        return out;
    }

    // allow for reduction operations
    template <class Function>
    constexpr Element reduce(Function f) const
        noexcept(noexcept(f(Element(), Element()))) {
        return f(c_[0],
                 f(c_[1],
                   f(c_[2],
                     f(c_[3], f(c_[4], f(c_[5], f(c_[6], f(c_[7], c_[8]))))))));
    }

    // ==========================================
    // Operators
    // ==========================================

    // equality comparison
    constexpr bool approx_eq(const Matrix& rhs) const noexcept {
        Element m = 0.0;
        Element d = 0.0;
        for (size_t i = 0; i < c_.size(); ++i) {
            m = std::max(std::abs(c_[i]), std::abs(rhs.c_[i]));
            d = std::abs(c_[i] - rhs.c_[i]);
            if (d > std::max(RELATIVE_TOL * m, ABSOLUTE_TOL)) {
                return false;
            }
        }
        return true;
    }

    constexpr bool operator==(const Matrix& rhs) const noexcept {
        return approx_eq(rhs);
    }
    constexpr bool operator!=(const Matrix& rhs) const noexcept {
        return !(this->operator==(rhs));
    }

    // element access
    Element& operator()(size_t i, size_t j) { return c_[3 * i + j]; }
    constexpr Element operator()(size_t i, size_t j) const {
        return c_[3 * i + j];
    }

    // unary +
    Matrix operator+() const { return *this; }

    // unary -
    Matrix operator-() const {
        Matrix out;
        for (size_t i = 0; i < c_.size(); ++i) {
            out.c_[i] = -c_[i];
        }
        return out;
    }

    GENERATE_MATRIX_OP(+)
    GENERATE_MATRIX_OP(-)
    GENERATE_MATRIX_OP(*)
    GENERATE_MATRIX_OP(/)

    template <class Stream>
    friend Stream& operator<<(Stream& stream, const Matrix& mat) {
        stream << "[";
        for (size_t i = 0; i < 3; ++i) {
            stream << (i == 0 ? "[" : " [");
            stream << mat(i, 0) << ", " << mat(i, 1) << ", " << mat(i, 2);
            stream << "]" << (i == 2 ? "" : "\n");
        }
        stream << "]";

        return stream;
    }

   private:
    std::array<Element, 9> c_;
};
