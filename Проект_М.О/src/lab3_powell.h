#pragma once
#include "oracle.h"
#include "svenn.h"
#include "lab2_search.h"
#include "lab3_nelder_mead.h"
#include <vector>
#include <cmath>
#include <algorithm>

struct PowellRow { int cycle; double a, b, f, shift, alpha1, alpha2; long n_f; };

struct PowellResult {
    double a, b, f;
    int cycles;
    long n_f;
    std::vector<PowellRow> rows;
};

// Лінійний пошук уздовж напрямку
inline void line_minimize(Oracle2& oracle, double& xa, double& xb,
                          double da, double db, double& f_x, double& alpha_best,
                          double delta = 0.1, double eps = 1e-4, double amax = 1000.0) {
    auto clampA = [](double v) { return std::min(std::max(v, 0.3), 5.0); };
    auto clampB = [](double v) { return std::min(std::max(v, 0.01), 1.0); };
    
    double ax = xa, bx = xb;
    Oracle line([&](double alpha) {
        return oracle.value(clampA(ax + alpha * da), clampB(bx + alpha * db));
    });
    
    auto bracket = svenn(line, 0.0, delta);
    auto best = golden_section(line, bracket.a, bracket.b, eps);
    
    alpha_best = best.x;
    xa = clampA(ax + best.x * da);
    xb = clampB(bx + best.x * db);
    f_x = oracle.value(xa, xb);
}

inline PowellResult powell(Oracle2& oracle, double a0, double b0,
                           double tol = 1e-4, int max_cycles = 50) {
    double xa = a0, xb = b0;
    double f_x = oracle.value(xa, xb);
    
    std::vector<std::pair<double, double>> dirs = {{1.0, 0.0}, {0.0, 1.0}};
    std::vector<PowellRow> rows;
    
    int cycle = 1;
    for (; cycle <= max_cycles; ++cycle) {
        double sa = xa, sb = xb, f_start = f_x;
        double gains[2], alphas[2];
        
        for (int i = 0; i < 2; ++i) {
            double current_f = f_x;
            line_minimize(oracle, xa, xb, dirs[i].first, dirs[i].second, f_x, alphas[i]);
            gains[i] = current_f - f_x;
        }
        
        double shift_a = xa - sa;
        double shift_b = xb - sb;
        double norm = std::sqrt(shift_a * shift_a + shift_b * shift_b);
        
        rows.push_back({cycle, xa, xb, f_x, norm, alphas[0], alphas[1], oracle.n_f});
        
        if (norm < tol || std::abs(f_start - f_x) <= 1e-12 * (1.0 + std::abs(f_x))) {
            break;
        }
        
        double dummy_alpha;
        line_minimize(oracle, xa, xb, shift_a / norm, shift_b / norm, f_x, dummy_alpha);
        
        int best_dir = (gains[0] > gains[1]) ? 0 : 1;
        dirs[best_dir] = {shift_a / norm, shift_b / norm};
        
        if (cycle % 3 == 0) {
            dirs = {{1.0, 0.0}, {0.0, 1.0}}; // Перезапуск напрямків
        }
    }
    return {xa, xb, f_x, cycle, oracle.n_f, rows};
}