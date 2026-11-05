#pragma once
#include "operators.h"
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

struct GradRow { int k; double E, grad, h; long n_f, n_grad, dist; };

struct GradResult {
    Image I;
    int iterations;
    long n_f, n_grad;
    double E, grad;
    bool converged;
    std::vector<GradRow> rows;
};

// kind: "gd" | "armijo" | "hb" | "nag"
inline GradResult gradient_method(const Image& noisy, double lam, double L,
                                  const std::string& kind, double beta = 0.0,
                                  double step = -1.0, int max_iter = 200,
                                  double tol = 1e-6, const Image* ref = nullptr) {
    double h = (step > 0.0) ? step : 1.0 / L;
    Image I = noisy, I_prev = noisy;
    long n_f = 0, n_grad = 0;
    double trial = h, g_norm = 0.0;
    bool converged = false;
    std::vector<GradRow> rows;
    
    for (int k = 1; k <= max_iter; ++k) {
        Image I_new(I.h, I.w);
        if (kind == "nag") {
            double b = (beta > 0.0) ? beta : double(k - 1) / (k + 2);
            Image Y(I.h, I.w);
            for (size_t t = 0; t < I.px.size(); ++t)
                Y.px[t] = I.px[t] + b * (I.px[t] - I_prev.px[t]);
            Image g = energy_grad(Y, noisy, lam); ++n_grad;
            g_norm = norm(g);
            for (size_t t = 0; t < I.px.size(); ++t) I_new.px[t] = Y.px[t] - h * g.px[t];
        } else {
            Image g = energy_grad(I, noisy, lam); ++n_grad;
            g_norm = norm(g);
            if (kind == "armijo") {
                double s = trial, f_I = energy(I, noisy, lam); ++n_f;
                double g2 = 0.0; for (double v : g.px) g2 += v * v;
                Image cand(I.h, I.w);
                for (int t = 0; t < 40; ++t) {
                    for (size_t q = 0; q < I.px.size(); ++q) cand.px[q] = I.px[q] - s * g.px[q];
                    double f_c = energy(cand, noisy, lam); ++n_f;
                    if (f_c <= f_I - 1e-4 * s * g2) break;
                    s *= 0.5;
                }
                I_new = cand; trial = s * 1.5; h = s;
            } else if (kind == "hb") {
                for (size_t t = 0; t < I.px.size(); ++t)
                    I_new.px[t] = I.px[t] - h * g.px[t] + beta * (I.px[t] - I_prev.px[t]);
            } else {
                for (size_t t = 0; t < I.px.size(); ++t) I_new.px[t] = I.px[t] - h * g.px[t];
            }
        }
        
        double e_k = energy(I, noisy, lam); ++n_f;
        long dist = (ref != nullptr) ? (long)(norm_diff(I, *ref) * 1e6) : 0; 
        rows.push_back({k - 1, e_k, g_norm, h, n_f, n_grad, dist});
        
        if (g_norm <= tol) { converged = true; break; }
        I_prev = I; I = I_new;
    }
    return {I, int(rows.size()) - 1, n_f, n_grad, energy(I, noisy, lam), g_norm, converged, rows};
}