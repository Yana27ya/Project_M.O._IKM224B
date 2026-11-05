#pragma once
#include "operators.h"
#include <vector>
#include <algorithm>
#include <iostream>

// Проєкція на відрізок [0, 1]
inline Image project(Image I, double lower = 0.0, double upper = 1.0) {
    for (double& v : I.px) v = std::min(std::max(v, lower), upper);
    return I;
}

struct Violation { double max_value; int count; };

inline Violation get_violation(const Image& I, double lower = 0.0, double upper = 1.0) {
    double worst = 0.0; int count = 0;
    for (double v : I.px) {
        double below = std::max(0.0, lower - v);
        double above = std::max(0.0, v - upper);
        worst = std::max(worst, std::max(below, above));
        if (v < lower || v > upper) ++count;
    }
    return {worst, count};
}

struct PGRow { int k; double E, G, viol; int active; long n_f, n_grad; };

struct PGResult {
    Image I; int iterations; long n_f, n_grad;
    double E, G; bool converged; std::vector<PGRow> rows;
};

// Метод проєкції градієнта
inline PGResult projected_gradient(const Image& raw, double lam, double L,
                                   bool accelerated = false,
                                   int max_iter = 200, double tol = 1e-6) {
    double h = 1.0 / L;
    Image I = project(raw), I_prev = I, Y = I;
    long n_f = 0, n_grad = 0;
    std::vector<PGRow> rows;
    bool converged = false;
    double g_map = 0.0;
    
    for (int k = 1; k <= max_iter; ++k) {
        const Image& point = accelerated ? Y : I;
        Image g = energy_grad(point, raw, lam); ++n_grad;
        Image I_new(I.h, I.w);
        
        for (size_t t = 0; t < I.px.size(); ++t) {
            I_new.px[t] = std::min(std::max(point.px[t] - h * g.px[t], 0.0), 1.0);
        }
        
        double acc = 0.0;
        for (size_t t = 0; t < I.px.size(); ++t) {
            double d = I_new.px[t] - point.px[t];
            acc += d * d;
        }
        g_map = std::sqrt(acc) / h; // градієнтне відображення
        
        int active = 0;
        for (double v : I.px) if (v <= 1e-12 || v >= 1.0 - 1e-12) ++active;
        Violation v = get_violation(I);
        double e_k = energy(I, raw, lam); ++n_f;
        
        rows.push_back({k - 1, e_k, g_map, v.max_value, active, n_f, n_grad});
        if (g_map <= tol) { converged = true; break; }
        
        if (accelerated) {
            double beta = double(k - 1) / (k + 2);
            for (size_t t = 0; t < I.px.size(); ++t)
                Y.px[t] = I_new.px[t] + beta * (I_new.px[t] - I_prev.px[t]);
            I_prev = I_new;
        }
        I = I_new;
    }
    return {I, int(rows.size()) - 1, n_f, n_grad, energy(I, raw, lam), g_map, converged, rows};
}


inline double penalty_value(const Image& I, double lower = 0.0, double upper = 1.0) {
    double s = 0.0;
    for (double v : I.px) {
        double below = std::max(0.0, lower - v), above = std::max(0.0, v - upper);
        s += below * below + above * above;
    }
    return s;
}

inline Image penalty_grad(const Image& I, double lower = 0.0, double upper = 1.0) {
    Image g(I.h, I.w);
    for (size_t t = 0; t < I.px.size(); ++t) {
        double v = I.px[t];
        g.px[t] = -2.0 * std::max(0.0, lower - v) + 2.0 * std::max(0.0, v - upper);
    }
    return g;
}

struct PenaltyRow { int s; double r, E, P, rP, violation; int outside; long n_f, n_grad; };
struct PenaltyResult { Image I; long n_f, n_grad; double E; std::vector<PenaltyRow> rows; };

inline PenaltyResult penalty_method(const Image& raw, double lam, double L, double r0 = 1.0, 
                                    double growth = 10.0, int outer = 7, int inner = 100, 
                                    double inner_tol = 1e-6, double tol = 1e-8) {
    Image I = raw;
    double r = r0;
    long n_f = 0, n_grad = 0;
    std::vector<PenaltyRow> rows;
    
    for (int s = 0; s < outer; ++s) {
        double L_r = L + 2.0 * r;
        double h = 1.0 / L_r;
        Image I_prev = I;
        
        for (int k = 1; k <= inner; ++k) {
            double beta = double(k - 1) / (k + 2);
            Image Y(I.h, I.w);
            for (size_t t = 0; t < I.px.size(); ++t) Y.px[t] = I.px[t] + beta * (I.px[t] - I_prev.px[t]);
            
            Image gE = energy_grad(Y, raw, lam);
            Image gP = penalty_grad(Y);
            ++n_grad;
            
            Image I_new(I.h, I.w);
            double g_norm_sq = 0.0;
            for (size_t t = 0; t < I.px.size(); ++t) {
                double total_g = gE.px[t] + r * gP.px[t];
                g_norm_sq += total_g * total_g;
                I_new.px[t] = Y.px[t] - h * total_g;
            }
            
            if (std::sqrt(g_norm_sq) <= inner_tol) {
                I_prev = I; I = I_new; break;
            }
            I_prev = I; I = I_new;
        }
        
        double P = penalty_value(I); ++n_f;
        Violation v = get_violation(I);
        rows.push_back({s, r, energy(I, raw, lam), P, r * P, v.max_value, v.count, n_f, n_grad});
        
        if (r * P <= tol) break;
        r *= growth;
    }
    return {I, n_f, n_grad, energy(I, raw, lam), rows};
}