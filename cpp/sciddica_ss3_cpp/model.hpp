#pragma once
#include <cstddef>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace sciddica {

template <typename T>
class Model {
public:
    static constexpr int ABI_VERSION = 1;
    Model(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {
        const std::size_t n = rows * cols;
        altitude_.assign(n, T(0));
        thickness_.assign(n, T(0));
        soil_depth_.assign(n, T(0));
        kinetic_head_.assign(n, T(0));
        momentum_x_.assign(n, T(0));
        momentum_y_.assign(n, T(0));
        barycentre_x_.assign(n, T(0));
        barycentre_y_.assign(n, T(0));
        max_impact_mask_.assign(n, (unsigned char)0);
    }
    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t size() const { return rows_ * cols_; }

    T* altitude() { return altitude_.data(); }
    T* thickness() { return thickness_.data(); }
    T* soil_depth() { return soil_depth_.data(); }
    T* kinetic_head() { return kinetic_head_.data(); }
    T* momentum_x() { return momentum_x_.data(); }
    T* momentum_y() { return momentum_y_.data(); }
    T* barycentre_x() { return barycentre_x_.data(); }
    T* barycentre_y() { return barycentre_y_.data(); }
    unsigned char* max_impact_mask() { return max_impact_mask_.data(); }

    const T* altitude() const { return altitude_.data(); }
    const T* thickness() const { return thickness_.data(); }

private:
    std::size_t rows_, cols_;
    std::vector<T> altitude_, thickness_, soil_depth_, kinetic_head_;
    std::vector<T> momentum_x_, momentum_y_, barycentre_x_, barycentre_y_;
    std::vector<unsigned char> max_impact_mask_;
};

using Model32 = Model<float>;
using Model64 = Model<double>;

} // namespace sciddica
