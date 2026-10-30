#include <iostream>
#include <cstdio>
#define NOMINMAX
#include <windows.h> 

#include "src/variant.h"
#include "src/image.h"
#include "src/rng.h"
#include "src/metrics.h"
#include "src/operators.h"
#include "src/lab4_gradient.h"
#include "src/lab5_constrained.h"

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    Image clean = clean_image(IMG_SIZE);
    Noisy result = add_noise(clean, NOISE_SIGMA, SEED);

Image raw_img = result.raw; 

    std::printf("=== ЧАСТИНА 5: УМОВНА ОПТИМІЗАЦІЯ ===\n");
    double lam = 0.624865;
    double L = 1.0 + 8.0 * lam;

    Violation v_raw = get_violation(raw_img);
    std::printf("Сире спостереження поза [0,1]: %d пікселів\n\n", v_raw.count);

    // 1. Проєкція градієнта
    std::printf("--- Проєкція градієнта ---\n");
    PGResult res_pg = projected_gradient(raw_img, lam, L);
    std::printf("Ітерацій: %d, Nf: %ld, Разом звернень: %ld\n", res_pg.iterations, res_pg.n_f, res_pg.n_f + res_pg.n_grad);
    std::printf("E* = %.8f, Порушення: %.8f\n\n", res_pg.E, res_pg.rows.back().viol);

    // 2. Метод штрафів
    std::printf("--- Метод зовнішніх штрафів ---\n");
    PenaltyResult res_pen = penalty_method(raw_img, lam, L);
    std::printf("Зовнішніх ітерацій: %d, Nf: %ld, Разом звернень: %ld\n", (int)res_pen.rows.size(), res_pen.n_f, res_pen.n_f + res_pen.n_grad);
    std::printf("E* = %.8f, Порушення: %.2e, Поза межами: %d пікс.\n\n", 
                res_pen.E, res_pen.rows.back().violation, res_pen.rows.back().outside);

    // Відновлення та вимір якості (з обмеженнями)
    write_pgm(res_pg.I, "lab5_restored.pgm");
    std::printf("PSNR з обмеженнями = %.4f дБ\n", psnr(res_pg.I, clean));

    return 0;
}