#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "rng.h"

struct Image {
    int h = 0, w = 0;
    std::vector<double> px;
    Image() = default;
    Image(int h_, int w_, double v = 0.0) : h(h_), w(w_), px(size_t(h_) * w_, v) {}
    double&       at(int i, int j)       { return px[size_t(i) * w + j]; }
    const double& at(int i, int j) const { return px[size_t(i) * w + j]; }
};

inline Image clean_image(int size = 128) {
    Image img(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double u = double(j) / (size - 1);
            double v = double(i) / (size - 1);
            double val = 0.25 + 0.35 * (0.6 * u + 0.4 * v);
            val += 0.08 * std::cos(2.0 * 3.14159265358979323846 * std::sqrt((u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5)));
            if ((u - 0.32) * (u - 0.32) + (v - 0.34) * (v - 0.34) <= 0.16 * 0.16) val = 0.92;
            if (std::fabs(u - 0.70) <= 0.15 && std::fabs(v - 0.30) <= 0.15) val = 0.10;
            if (v > 0.60 && v < 0.90 && std::fabs(u - 0.30) < 0.9 * (v - 0.60)) val = 0.68;
            if (u > 0.58 && u < 0.92 && v > 0.60 && v < 0.90) {
                int cell = std::max(3, size / 12);
                val = (((j / cell) + (i / cell)) % 2 == 0) ? 0.85 : 0.20;
            }
            img.at(i, j) = std::min(std::max(val, 0.0), 1.0);
        }
    }
    for (int idx = 0; idx < 4; ++idx) {
        int thickness = std::max(2, (4 - idx) * std::max(1, size / 64));
        int row = int(0.04 * size) + idx * std::max(3, size / 22);
        if (row + thickness < size)
            for (int i = row; i < row + thickness; ++i)
                for (int j = int(0.05 * size); j < int(0.95 * size); ++j)
                    img.at(i, j) = 0.98;
    }
    return img;
}

inline void write_pgm(const Image& img, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("не вдалося створити " + path);
    out << "P5\n" << img.w << " " << img.h << "\n255\n";
    for (double v : img.px) {
        double c = std::min(std::max(v, 0.0), 1.0);
        out.put(static_cast<char>(int(std::lround(c * 255.0))));
    }
}

struct Noisy {
    Image raw;
    Image noisy;
};

inline Noisy add_noise(const Image& clean, double sigma, uint64_t seed) {
    Rng rng(seed);
    Image raw(clean.h, clean.w), noisy(clean.h, clean.w);
    for (int i = 0; i < clean.h; ++i) {
        for (int j = 0; j < clean.w; ++j) {
            double v = clean.at(i, j) + sigma * rng.normal();
            raw.at(i, j)   = v;
            noisy.at(i, j) = std::min(std::max(v, 0.0), 1.0);
        }
    }
    return {raw, noisy};
}