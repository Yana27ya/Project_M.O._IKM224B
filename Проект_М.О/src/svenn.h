#pragma once
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "oracle.h"

struct SvennInterval {
    double a;
    double b;
};

inline SvennInterval svenn(Oracle& oracle, double x0, double delta) {
    std::printf("k\tx_k\t\tf(x_k)\t\tdelta\n");
    std::printf("--------------------------------------------------\n");

    double x_k = x0;
    double f_k = oracle.value(x_k);
    
    double f_plus = oracle.value(x_k + delta);
    double f_minus = oracle.value(x_k - delta);
    
    double step = 0;
    
    if (f_minus >= f_k && f_k >= f_plus) {
        step = delta;
    } else if (f_minus <= f_k && f_k <= f_plus) {
        step = -delta;
    } else if (f_minus >= f_k && f_k <= f_plus) {
        std::printf("0\t%.4f\t\t%.6f\t%.4f\n", x_k, f_k, delta);
        return {x_k - std::abs(delta), x_k + std::abs(delta)};
    } else {
        std::printf("Помилка: функція не є унімодальною або крок завеликий.\n");
        return {x_k, x_k};
    }

    int k = 0;
    double x_prev = x_k;
    double f_prev = f_k;
    
    x_k = x0 + step;
    f_k = (step > 0) ? f_plus : f_minus;
    
    std::printf("%d\t%.4f\t\t%.6f\t%.4f\n", k, x_prev, f_prev, step);
    
    while (true) {
        k++;
        step *= 2.0;
        double x_next = x_k + step;
        double f_next = oracle.value(x_next);
        
        std::printf("%d\t%.4f\t\t%.6f\t%.4f\n", k, x_k, f_k, step);
        
        if (f_next >= f_k) {
            std::printf("%d\t%.4f\t\t%.6f\t-\n", k+1, x_next, f_next);
            double a = std::min(x_prev, x_next);
            double b = std::max(x_prev, x_next);
            return {a, b};
        }
        
        x_prev = x_k;
        f_prev = f_k;
        x_k = x_next;
        f_k = f_next;
    }
}