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
inline double total_variation(const Image& I) {
    double tv = 0.0;
    for (int y = 0; y < I.h - 1; ++y) {
        for (int x = 0; x < I.w - 1; ++x) {
            double dx = I.px[y * I.w + (x + 1)] - I.px[y * I.w + x];
            double dy = I.px[(y + 1) * I.w + x] - I.px[y * I.w + x];
            tv += std::sqrt(dx * dx + dy * dy);
        }
    }
    return tv;
}


inline double ssim(const Image& I1, const Image& I2) {
    if (I1.px.size() != I2.px.size()) return 0.0;
    
    int n = I1.px.size();
    double mu1 = 0.0, mu2 = 0.0;
    
    // Рахуємо середнє значення
    for (int i = 0; i < n; ++i) {
        mu1 += I1.px[i];
        mu2 += I2.px[i];
    }
    mu1 /= n;
    mu2 /= n;

    double var1 = 0.0, var2 = 0.0, covar = 0.0;
    // Рахуємо дисперсію та коваріацію
    for (int i = 0; i < n; ++i) {
        var1 += (I1.px[i] - mu1) * (I1.px[i] - mu1);
        var2 += (I2.px[i] - mu2) * (I2.px[i] - mu2);
        covar += (I1.px[i] - mu1) * (I2.px[i] - mu2);
    }
    var1 /= n;
    var2 /= n;
    covar /= n;

    double c1 = 0.0001; 
    double c2 = 0.0009; 

    double numerator = (2 * mu1 * mu2 + c1) * (2 * covar + c2);
    double denominator = (mu1 * mu1 + mu2 * mu2 + c1) * (var1 + var2 + c2);
    
    return numerator / denominator;
}