#pragma once
#include <cstdint>

constexpr int      STUDENT_ID   = 7;
constexpr uint64_t SEED         = (uint64_t(STUDENT_ID) * 7919 + 104729) % 2147483647ULL;
constexpr double   NOISE_SIGMA  = 0.1743;
constexpr double   X0           = 0.5663;
constexpr double   DELTA0       = 0.2026;
constexpr double   EPS          = 1e-4;
constexpr double   LAMBDA0      = 0.6319;
constexpr double   SIMPLEX_EDGE = 0.6280;
constexpr int      IMG_SIZE     = 128; 