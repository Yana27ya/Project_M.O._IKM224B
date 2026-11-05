#pragma once
#include <cstdint>
#include <cmath>

struct Rng {
    uint64_t state;
    explicit Rng(uint64_t seed) : state(seed) {}

    uint64_t next_u64() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state;
    }
    double uniform() {
        return (next_u64() >> 11) * (1.0 / 9007199254740992.0);
    }
    double normal() {
        double u1 = uniform();
        if (u1 < 1e-300) u1 = 1e-300;
        double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
    }
};