#pragma once
#include "oracle.h"
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

const double TAU = (std::sqrt(5.0) - 1.0) / 2.0; 

struct SearchRow {
    int k; double a, b, x1, x2, f1, f2, L, new_point; std::string decision;
};

struct SearchResult {
    double x, f, a, b;
    int iterations;
    long n_f;
    std::vector<SearchRow> rows;
};

inline SearchResult dichotomy(Oracle& oracle, double a, double b,
                              double eps = 1e-5, double delta = -1.0,
                              int max_iter = 500) {
    if (delta < 0.0) delta = eps / 3.0;
    if (!(delta > 0.0 && delta < eps))
        throw std::invalid_argument("має виконуватися 0 < delta < eps");
    std::vector<SearchRow> rows;
    int k = 0;
    while ((b - a) > eps && k < max_iter) {
        double x1 = (a + b - delta) / 2.0;
        double x2 = (a + b + delta) / 2.0;
        double f1 = oracle.value(x1);                  // звернення 1
        double f2 = oracle.value(x2);                  // звернення 2
        std::string decision;
        if (f1 <= f2) { decision = "f(x1) <= f(x2) -> b = x2"; b = x2; }
        else          { decision = "f(x1) >  f(x2) -> a = x1"; a = x1; }
        ++k;
        rows.push_back({k, a, b, x1, x2, f1, f2, b - a, 0.0, decision});
    }
    double x_best = (a + b) / 2.0;
    return {x_best, oracle.value(x_best), a, b, k, oracle.n_f, rows};
}

inline SearchResult golden_section(Oracle& oracle, double a, double b,
                                   double eps = 1e-5, int max_iter = 500) {
    double x1 = b - TAU * (b - a);
    double x2 = a + TAU * (b - a);
    double f1 = oracle.value(x1);                      // два звернення лише
    double f2 = oracle.value(x2);                      // на самому початку
    std::vector<SearchRow> rows;
    int k = 0;
    while ((b - a) > eps && k < max_iter) {
        std::string decision;
        double new_point;
        if (f1 <= f2) {
            decision = "f(x1) <= f(x2) -> b = x2, успадкували x1";
            b = x2;
            x2 = x1; f2 = f1;                          // <- фокус
            x1 = b - TAU * (b - a);
            f1 = oracle.value(x1);                     // єдине нове звернення
            new_point = x1;
        } else {
            decision = "f(x1) >  f(x2) -> a = x1, успадкували x2";
            a = x1;
            x1 = x2; f1 = f2;                          // <- фокус
            x2 = a + TAU * (b - a);
            f2 = oracle.value(x2);                     // єдине нове звернення
            new_point = x2;
        }
        ++k;
        rows.push_back({k, a, b, x1, x2, f1, f2, b - a, new_point, decision});
    }
    double x_best = (a + b) / 2.0;
    return {x_best, oracle.value(x_best), a, b, k, oracle.n_f, rows};
}