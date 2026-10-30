#pragma once
#include "image.h"
#include <cmath>
#include <limits>

inline double mse(const Image& a, const Image& b) {
    double s = 0.0;
    for (size_t k = 0; k < a.px.size(); ++k) {
        double d = a.px[k] - b.px[k];
        s += d * d;
    }
    return s / double(a.px.size());
}

inline double psnr(const Image& a, const Image& b) {
    double v = mse(a, b);
    if (v <= 0.0) return std::numeric_limits<double>::infinity();
    return 10.0 * std::log10(1.0 / v);
}