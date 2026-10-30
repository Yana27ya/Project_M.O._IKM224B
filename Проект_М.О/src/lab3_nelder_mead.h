#pragma once
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <functional>

struct NMRow { int k; std::string op; double a, b, f, size; long n_f; };

struct NMResult {
    double a, b, f;
    int iterations;
    long n_f;
    std::vector<NMRow> rows;
};

// Оракул для функції двох змінних
struct Oracle2 {
    std::function<double(double, double)> func;
    long n_f = 0;
    explicit Oracle2(std::function<double(double, double)> f) : func(std::move(f)) {}
    double value(double a, double b) { ++n_f; return func(a, b); }
    double peek(double a, double b) const { return func(a, b); }
    void reset() { n_f = 0; }
};

inline NMResult nelder_mead(Oracle2& oracle, double a0, double b0, double step,
                            double tol = 1e-7, int max_iter = 200) {
    const double ALPHA = 1.0, GAMMA = 2.0, RHO = 0.5, SIGMA = 0.5;
    auto clampA = [](double v) { return std::min(std::max(v, 0.3), 5.0); };
    auto clampB = [](double v) { return std::min(std::max(v, 0.01), 1.0); };
    
    std::array<double, 3> xa{a0, clampA(a0 + step), a0};
    std::array<double, 3> xb{b0, b0, clampB(b0 + 0.4 * step)};
    std::array<double, 3> fv{};
    
    for (int i = 0; i < 3; ++i) fv[size_t(i)] = oracle.value(xa[size_t(i)], xb[size_t(i)]);
    std::vector<NMRow> rows;
    int k = 1;
    
    for (; k <= max_iter; ++k) {
        std::array<int, 3> idx{0, 1, 2};
        std::sort(idx.begin(), idx.end(), [&](int p, int q) { return fv[size_t(p)] < fv[size_t(q)]; });
        std::array<double, 3> na{}, nb{}, nf{};
        for (int i = 0; i < 3; ++i) {
            na[size_t(i)] = xa[size_t(idx[size_t(i)])];
            nb[size_t(i)] = xb[size_t(idx[size_t(i)])];
            nf[size_t(i)] = fv[size_t(idx[size_t(i)])];
        }
        xa = na; xb = nb; fv = nf;
        
        double ca = (xa[0] + xa[1]) / 2.0, cb = (xb[0] + xb[1]) / 2.0;
        double ra = clampA(ca + ALPHA * (ca - xa[2]));
        double rb = clampB(cb + ALPHA * (cb - xb[2]));
        double fr = oracle.value(ra, rb);
        std::string op;
        
        if (fr < fv[0]) {
            double ea = clampA(ca + GAMMA * (ra - ca)), eb = clampB(cb + GAMMA * (rb - cb));
            double fe = oracle.value(ea, eb);
            if (fe < fr) { xa[2] = ea; xb[2] = eb; fv[2] = fe; op = "розтягування"; }
            else         { xa[2] = ra; xb[2] = rb; fv[2] = fr; op = "відображення"; }
        } else if (fr < fv[1]) {
            xa[2] = ra; xb[2] = rb; fv[2] = fr; op = "відображення";
        } else {
            double ka, kb;
            if (fr < fv[2]) { ka = ca + RHO * (ra - ca); kb = cb + RHO * (rb - cb); op = "зовнішнє стиснення"; }
            else            { ka = ca + RHO * (xa[2] - ca); kb = cb + RHO * (xb[2] - cb); op = "внутрішнє стиснення"; }
            ka = clampA(ka); kb = clampB(kb);
            double fc = oracle.value(ka, kb);
            
            if (fc < std::min(fr, fv[2])) { xa[2] = ka; xb[2] = kb; fv[2] = fc; }
            else {
                for (int i = 1; i < 3; ++i) {
                    xa[size_t(i)] = clampA(xa[0] + SIGMA * (xa[size_t(i)] - xa[0]));
                    xb[size_t(i)] = clampB(xb[0] + SIGMA * (xb[size_t(i)] - xb[0]));
                    fv[size_t(i)] = oracle.value(xa[size_t(i)], xb[size_t(i)]);
                }
                op = "редукція";
            }
        }
        
        int ib = int(std::min_element(fv.begin(), fv.end()) - fv.begin());
        double size = 0.0;
        for (int i = 0; i < 3; ++i) {
            double da = xa[size_t(i)] - xa[size_t(ib)], db = xb[size_t(i)] - xb[size_t(ib)];
            size = std::max(size, std::sqrt(da * da + db * db));
        }
        rows.push_back({k, op, xa[size_t(ib)], xb[size_t(ib)], fv[size_t(ib)], size, oracle.n_f});
        double spread = *std::max_element(fv.begin(), fv.end()) - *std::min_element(fv.begin(), fv.end());
        if (size <= tol && spread <= tol) break;
    }
    int ib = int(std::min_element(fv.begin(), fv.end()) - fv.begin());
    return {xa[size_t(ib)], xb[size_t(ib)], fv[size_t(ib)], k, oracle.n_f, rows};
}