#pragma once
#include "image.h"
#include <cmath>
#include <vector>

struct Field { Image x, y; }; 
inline Field grad_xy(const Image& u) {
    Field g{Image(u.h, u.w), Image(u.h, u.w)};
    for (int i = 0; i < u.h; ++i)
        for (int j = 0; j < u.w; ++j) {
            g.x.at(i, j) = (j < u.w - 1) ? u.at(i, j + 1) - u.at(i, j) : 0.0;
            g.y.at(i, j) = (i < u.h - 1) ? u.at(i + 1, j) - u.at(i, j) : 0.0;
        }
    return g;
}

inline Image divergence(const Image& px, const Image& py) {
    Image d(px.h, px.w);
    for (int i = 0; i < px.h; ++i)
        for (int j = 0; j < px.w; ++j) {
            double value = 0.0;
            if (j == 0)                 value += px.at(i, 0);
            else if (j < px.w - 1)      value += px.at(i, j) - px.at(i, j - 1);
            else                        value -= px.at(i, px.w - 2);
            
            if (i == 0)                 value += py.at(0, j);
            else if (i < px.h - 1)      value += py.at(i, j) - py.at(i - 1, j);
            else                        value -= py.at(px.h - 2, j);
            d.at(i, j) = value;
        }
    return d;
}

inline Image laplacian(const Image& u) {
    Field g = grad_xy(u);
    return divergence(g.x, g.y);
}

inline double energy(const Image& I, const Image& noisy, double lam) {
    Field g = grad_xy(I);
    double fit = 0.0, smooth = 0.0;
    for (size_t k = 0; k < I.px.size(); ++k) {
        double d = I.px[k] - noisy.px[k];
        fit += d * d;
        smooth += g.x.px[k] * g.x.px[k] + g.y.px[k] * g.y.px[k];
    }
    return 0.5 * fit + 0.5 * lam * smooth;
}

inline Image energy_grad(const Image& I, const Image& noisy, double lam) {
    Image lap = laplacian(I);
    Image g(I.h, I.w);
    for (size_t k = 0; k < I.px.size(); ++k)
        g.px[k] = (I.px[k] - noisy.px[k]) - lam * lap.px[k];
    return g;
}

inline double norm(const Image& a) {
    double s = 0.0;
    for (double v : a.px) s += v * v;
    return std::sqrt(s);
}

// Норма різниці між двома зображеннями
inline double norm_diff(const Image& a, const Image& b) {
    double s = 0.0;
    for (size_t k = 0; k < a.px.size(); ++k) {
        double d = a.px[k] - b.px[k];
        s += d * d;
    }
    return std::sqrt(s);
}