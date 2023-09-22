#include <cmath>
#include <vector>
#include <memory>
#include <cassert>
#include <fstream>
#include <cstring>
#include <iostream>


template <typename T>
class CubicHermiteSpline {
public:
    CubicHermiteSpline(const T * x_ptr,
                       const T * y_ptr,
                       const T * m_ptr,
                       const size_t size):
        x_(x_ptr, x_ptr + size),
        y_(y_ptr, y_ptr + size),
        m_(m_ptr, m_ptr + size),
        size_(size) {};

    T get_interpolated_value(const T x) const; 

private:
    bool is_x_in_boundary(const size_t idx, const T x) const; 

    size_t binary_search_(const T x) const; 

    T interp_func_(const T t, const size_t idx) const; 

    std::vector<T> x_, y_, m_;
    size_t size_;
};


template <typename T>
class MonotoneCubicInterpolation {
public:
    MonotoneCubicInterpolation(const T * x_ptr, const T * y_ptr, const size_t size) {}; 

    T operator()(const T x) const {}; 

private:
    static constexpr const T keps = 1e-10;
    std::unique_ptr<CubicHermiteSpline<T>> spliner_ptr_;
};