#pragma once

#include <unordered_set>
#include <vector>

class NeighborList {
   public:
    using Element = Set<size_t>;
    using Container = List<Element>;

    NeighborList(size_t size) : c_(size) {}

    void resize(size_t size) { c_.resize(size); }

    // Sets 'first' and 'second' to be neighbors, symmetrically
    void set(size_t first, size_t second) noexcept {
        if (first != second) {
            c_[first].insert(second);
            c_[second].insert(first);
        }
    }

    void erase(size_t first, size_t second) noexcept {
        if (first != second) {
            c_[first].erase(second);
            c_[second].erase(first);
        }
    }

    // returns set of neighbors of argument
    const Element& operator[](size_t id) noexcept {
        return c_[id];
    }

    // returns true if first is a neighbor of second
    // Can only be used for queries! Use 'set' for setting neighbors
    bool operator()(size_t first, size_t second) const {
        return c_[first].count(second) > 0;
    }

   private:
    Container c_;
};
