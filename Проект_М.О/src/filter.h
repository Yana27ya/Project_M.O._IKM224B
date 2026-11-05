#pragma once
#include "image.h"
#include <cmath>
#include <vector>
#include <algorithm>

inline Image apply_gaussian_blur(const Image& img, double sigma) {
    if (sigma <= 0.01) return img; 

    int radius = std::ceil(3.0 * sigma);
    int size = 2 * radius + 1;
    std::vector<double> kernel(size);
    double sum = 0.0;

    // створюємо ядро Гаусса
    for (int i = -radius; i <= radius; ++i) {
        kernel[i + radius] = std::exp(-(i * i) / (2.0 * sigma * sigma));
        sum += kernel[i + radius];
    }
    for (int i = 0; i < size; ++i) kernel[i] /= sum;

    Image temp(img.h, img.w);
    // розмиваємо по горизонталі
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            double val = 0.0;
            for (int k = -radius; k <= radius; ++k) {
                int col = std::min(std::max(j + k, 0), img.w - 1);
                val += img.at(i, col) * kernel[k + radius];
            }
            temp.at(i, j) = val;
        }
    }

    Image result(img.h, img.w);
    // розмиваємо по вертикалі
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            double val = 0.0;
            for (int k = -radius; k <= radius; ++k) {
                int row = std::min(std::max(i + k, 0), img.h - 1);
                val += temp.at(row, j) * kernel[k + radius];
            }
            result.at(i, j) = val;
        }
    } 
    return result;
}

    // Дзеркальне відбиття індексу
inline int mirror(int index, int size) {
    if (index < 0) index = -index - 1;
    if (index >= size) index = 2 * size - index - 1;
    return std::min(std::max(index, 0), size - 1);
}

// Білатеральний фільтр
inline Image bilateral_filter(const Image& img, double sigma_spatial, double sigma_range) {
    sigma_spatial = std::max(sigma_spatial, 1e-3);
    sigma_range = std::max(sigma_range, 1e-3);
    int r = int(std::min(std::max(std::lround(2.0 * sigma_spatial), 1L), 5L));
    
    Image out(img.h, img.w);
    double two_ss = 2.0 * sigma_spatial * sigma_spatial;
    double two_sr = 2.0 * sigma_range * sigma_range;
    
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            double center = img.at(i, j);
            double num = 0.0, den = 0.0;
            
            for (int dy = -r; dy <= r; ++dy) {
                for (int dx = -r; dx <= r; ++dx) {
                    double neighbour = img.at(mirror(i + dy, img.h), mirror(j + dx, img.w));
                    double ws = std::exp(-(double(dx) * dx + double(dy) * dy) / two_ss);
                    double diff = neighbour - center;
                    double wr = std::exp(-(diff * diff) / two_sr);
                    double weight = ws * wr;
                    
                    num += weight * neighbour;
                    den += weight;
                }
            }
            out.at(i, j) = num / den;
        }
    }
    return out;
}