#pragma once
#include <functional>
#include <vector>
#include <utility>

struct Oracle {
    std::function<double(double)> func;
    long n_f = 0;
    std::vector<std::pair<double, double>> history;

    explicit Oracle(std::function<double(double)> f) : func(std::move(f)) {}

    double value(double x) {
        ++n_f;
        double v = func(x);
        history.emplace_back(x, v);
        return v;
    }
    double peek(double x) const { return func(x); }
    void reset() { n_f = 0; history.clear(); }
};