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

#define GENERATE_VECTOR_OP(OP)                                           \
    constexpr Vector& operator OP##=(const Vector& rhs) {                \
        c_[0] OP## = rhs.c_[0];                                          \
        c_[1] OP## = rhs.c_[1];                                          \
        c_[2] OP## = rhs.c_[2];                                          \
        return *this;                                                    \
    }                                                                    \
    constexpr Vector& operator OP##=(Element rhs) {                      \
        c_[0] OP## = rhs;                                                \
        c_[1] OP## = rhs;                                                \
        c_[2] OP## = rhs;                                                \
        return *this;                                                    \
    }                                                                    \
    friend constexpr Vector operator OP(Vector lhs, const Vector& rhs) { \
        lhs OP## = rhs;                                                  \
        return lhs;                                                      \
    }                                                                    \
    friend constexpr Vector operator OP(Vector lhs, Element rhs) {       \
        lhs OP## = rhs;                                                  \
        return lhs;                                                      \
    }                                                                    \
    friend constexpr Vector operator OP(Element lhs, Vector rhs) {       \
        rhs.c_[0] = lhs OP rhs.c_[0];                                    \
        rhs.c_[1] = lhs OP rhs.c_[1];                                    \
        rhs.c_[2] = lhs OP rhs.c_[2];                                    \
        return rhs;                                                      \
    }

class Vector {
   public:
    using Element = Scalar;
    using Container = Array<Element, 3>;

    // STL-style aliases
    using value_type = typename Container::value_type;
    using reference = typename Container::reference;
    using const_reference = typename Container::const_reference;
    using pointer = typename Container::pointer;
    using const_pointer = typename Container::const_pointer;
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;

    // explicitly default normal constructors/assignment
    constexpr Vector() noexcept : c_{} {}
    Vector(const Vector&) = default;
    Vector(Vector&&) = default;
    Vector& operator=(const Vector&) = default;
    Vector& operator=(Vector&&) = default;

    // construct from list
    constexpr Vector(Element x, Element y, Element z) noexcept : c_{x, y, z} {}

    // ==========================================
    // Iterator functions
    // ==========================================

    iterator begin() { return c_.begin(); }
    const_iterator begin() const { return c_.begin(); }
    iterator end() { return c_.end(); }
    const_iterator end() const { return c_.end(); }

    // ==========================================
    // Mathematical functions
    // ==========================================

    // Compute the squared-Euclidean length of the vector
    constexpr Element norm2() const {
        return c_[0] * c_[0] + c_[1] * c_[1] + c_[2] * c_[2];
    }

    // Compute the L2-norm (Euclidean length) of the vector
    // can't make constexpr because std::sqrt isn't constexpr :(
    Element norm() const noexcept { return std::sqrt(norm2()); }

    // normalize vector in-place, and return the computed magnitude of the
    // vector pre-normalization.
    Element normalize() {
        Element n = norm();
        // make sure not to normalize 0 vectors!
        if (n > 0.0) {
            c_[0] /= n;
            c_[1] /= n;
            c_[2] /= n;
        }
        return n;
    }

    // allow for elementwise mapping of a function
    template <class Function>
    constexpr Vector map(Function f) const noexcept(noexcept(f(Element()))) {
        Vector out;
        for (size_t i = 0; i < c_.size(); ++i) {
            out(i) = f(c_[i]);
        }
        return out;
    }

    // allow for reduction/fold operations
    template <class Function>
    constexpr Element reduce(Function f) const
        noexcept(noexcept(f(Element(), Element()))) {
        return f(c_[0], f(c_[1], c_[2]));
    }

    // ==========================================
    // Operators
    // ==========================================

    // element access
    inline Element& operator()(size_t i) { return c_[i]; }
    constexpr Element operator()(size_t i) const { return c_[i]; }
    inline Element& operator[](size_t i) { return c_[i]; }
    constexpr Element operator[](size_t i) const { return c_[i]; }

    inline Element& x() { return c_[0]; }
    inline Element& y() { return c_[1]; }
    inline Element& z() { return c_[2]; }
    constexpr Element x() const { return c_[0]; }
    constexpr Element y() const { return c_[1]; }
    constexpr Element z() const { return c_[2]; }

    // equality comparison
    constexpr bool approx_eq(const Vector& rhs) const noexcept {
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

    constexpr bool operator==(const Vector& rhs) const noexcept {
        return approx_eq(rhs);
    }
    constexpr bool operator!=(const Vector& rhs) const noexcept {
        return !(this->operator==(rhs));
    }

    // unary +
    constexpr Vector operator+() const { return *this; }

    // unary -
    constexpr Vector operator-() const {
        Vector out;
        out.c_[0] = -c_[0];
        out.c_[1] = -c_[1];
        out.c_[2] = -c_[2];
        return out;
    }

    // all arithmetic operators
    GENERATE_VECTOR_OP(+)
    GENERATE_VECTOR_OP(-)
    GENERATE_VECTOR_OP(*)
    GENERATE_VECTOR_OP(/)

    // specialization of '%' operator for cross products
    friend constexpr Vector operator%(const Vector& lhs, const Vector& rhs) {
        return {lhs.c_[1] * rhs.c_[2] - lhs.c_[2] * rhs.c_[1],
                lhs.c_[2] * rhs.c_[0] - lhs.c_[0] * rhs.c_[2],
                lhs.c_[0] * rhs.c_[1] - lhs.c_[1] * rhs.c_[0]};
    }

    // allow for human-readable output
    template <class Stream>
    friend inline Stream& operator<<(Stream& stream, const Vector& v) {
        stream << "[" << v.c_[0] << ", " << v.c_[1] << ", " << v.c_[2] << "]";
        return stream;
    }

   private:
    Container c_;
};

#undef GENERATE_VECTOR_OP
