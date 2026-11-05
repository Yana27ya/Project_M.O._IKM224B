#include <iostream>
#include <cstdio>
#define NOMINMAX
#include <windows.h> 
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "src/variant.h"
#include "src/image.h"
#include "src/rng.h"
#include "src/metrics.h"
#include "src/operators.h"
#include "src/lab4_gradient.h"
#include "src/lab2_search.h"
#include "src/oracle.h"
#include "src/svenn.h"
#include "src/filter.h"
#include "src/lab3_nelder_mead.h"
#include "src/lab3_powell.h"
#include "src/lab5_constrained.h"

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    std::printf("ЧАСТИНА 0: ГЕНЕРАЦІЯ ЗОБРАЖЕННЯ \n");
    Image clean = clean_image(IMG_SIZE);
    Noisy result = add_noise(clean, NOISE_SIGMA, SEED);
    Image raw_img = result.raw; 
    
    write_pgm(clean, "original.pgm");
    write_pgm(result.noisy, "noisy.pgm");
    std::printf("PSNR зашумленого зображення: %.4f дБ\n\n", psnr(result.noisy, clean));

    double sigma_star = 0.0; 

    // лаб 1 та 2
    {
        auto blur_mse_func = [&](double sigma) {
            return mse(apply_gaussian_blur(result.noisy, sigma), clean);
        };
        Oracle oracle(blur_mse_func);

        std::printf(" ЛАБОРАТОРНА 1: АЛГОРИТМ СВЕННА \n");
        std::printf("Початкова точка x0 = %.4f, крок delta0 = %.4f\n", X0, DELTA0);
        SvennInterval bracket = svenn(oracle, X0, DELTA0);
        std::printf("Відрізок локалізації: [%.4f, %.4f]\n", bracket.a, bracket.b);
        std::printf("Звернень до оракула: %ld\n\n", oracle.n_f);

        std::printf(" ЛАБОРАТОРНА 2: ДИХОТОМІЯ ТА ЗОЛОТИЙ ПЕРЕРІЗ \n");
        oracle.reset(); 
        SearchResult res_d = dichotomy(oracle, bracket.a, bracket.b, EPS);
        
        oracle.reset(); 
        SearchResult res_g = golden_section(oracle, bracket.a, bracket.b, EPS);

        std::printf("Дихотомія        sigma* = %.6f  ітерацій = %2d  N_f = %ld\n", res_d.x, res_d.iterations, res_d.n_f);
        std::printf("Золотий переріз  sigma* = %.6f  ітерацій = %2d  N_f = %ld\n", res_g.x, res_g.iterations, res_g.n_f);

        sigma_star = res_g.x; 
        Image restored = apply_gaussian_blur(result.noisy, sigma_star);
        write_pgm(restored, "lab2_restored.pgm"); 

        std::printf("\n-Відновлення зображення (Фільтр Гаусса) \n");
        std::printf("було:  MSE = %.8f  PSNR = %.4f дБ\n", mse(result.noisy, clean), psnr(result.noisy, clean));
        std::printf("стало: MSE = %.8f  PSNR = %.4f дБ\n\n", mse(restored, clean), psnr(restored, clean));
    }

    // лаб 3
    {
        std::printf(" ЛАБОРАТОРНА 3: НЕЛДЕР-МІД ТА ПАУЕЛЛ \n");
        double a0 = std::min(std::max(2.0 * sigma_star, 0.5), 4.0);
        double b0 = std::min(std::max(0.15 + 0.1 * NOISE_SIGMA, 0.02), 0.9);
        std::printf("Початкова точка theta0 = (%.4f, %.4f)\n", a0, b0);
        
        auto objective2 = [&](double a, double b) {
            auto clampA = [](double v) { return std::min(std::max(v, 0.3), 5.0); };
            auto clampB = [](double v) { return std::min(std::max(v, 0.01), 1.0); };
            double ca = clampA(a), cb = clampB(b);
            double val = mse(bilateral_filter(result.noisy, ca, cb), clean);
            double da = a - ca, db = b - cb;
            return val + 1000.0 * (da * da + db * db);
        };
        
        Oracle2 oracle2(objective2);
        std::printf("f(theta0) = %.8f\n\n", oracle2.peek(a0, b0));

        oracle2.reset();
        NMResult nm_res = nelder_mead(oracle2, a0, b0, 0.628, 1e-12);
        
        oracle2.reset();
        PowellResult pw_res = powell(oracle2, a0, b0);
        
        std::printf("Нелдер-Мід: theta* = (%.6f, %.6f)  ітер = %2d  N_f = %ld\n", nm_res.a, nm_res.b, nm_res.iterations, nm_res.n_f);
        std::printf("Пауелл:     theta* = (%.6f, %.6f)  циклів = %2d  N_f = %ld\n", pw_res.a, pw_res.b, pw_res.cycles, pw_res.n_f);
        
        double best_a = nm_res.f < pw_res.f ? nm_res.a : pw_res.a;
        double best_b = nm_res.f < pw_res.f ? nm_res.b : pw_res.b;
        Image restored = bilateral_filter(result.noisy, best_a, best_b);
        write_pgm(restored, "lab3_restored.pgm");
        
        std::printf("\n Відновлення зображення (Білатеральний фільтр) \n");
        std::printf("PSNR = %.4f дБ\n\n", psnr(restored, clean));

        std::ofstream f_grid("lab3_grid.txt");
        for (int i = 0; i <= 50; ++i) {
            double ss = 0.5 + i * (4.5 / 50.0);
            for (int j = 0; j <= 50; ++j) {
                double sr = 0.01 + j * (0.99 / 50.0);
                f_grid << ss << " " << sr << " " << oracle2.peek(ss, sr) << "\n";
            }
        }
        f_grid.close();
        
        std::ofstream f_nm("lab3_nm.txt");
        f_nm << a0 << " " << b0 << " " << nm_res.rows[0].f << "\n"; 
        for (const auto& row : nm_res.rows) f_nm << row.a << " " << row.b << " " << row.f << "\n";
        f_nm.close();
        
        std::ofstream f_pw("lab3_powell.txt");
        f_pw << a0 << " " << b0 << " " << pw_res.rows[0].f << "\n"; 
        for (const auto& row : pw_res.rows) f_pw << row.a << " " << row.b << " " << row.f << "\n";
        f_pw.close();
    }

    // лаб 4
    double lam = 0.624865;
    double L = 6.0;

    {
        std::printf(" ЧАСТИНА 4: ГРАДІЄНТНІ МЕТОДИ \n");

        std::printf("Шукаємо оптимальне lambda* (Свенн + Золотий переріз)...\n");
        Oracle oracle_lam([&](double l) {
            Image temp = gradient_method(result.noisy, l, 1.0 + 8.0 * l, "gd", 0.0, -1.0, 50, 1e-4).I;
            return mse(temp, clean);
        });
        
        auto svenn_res = svenn(oracle_lam, 0.1, 0.1); 
        auto gs_res = golden_section(oracle_lam, svenn_res.a, svenn_res.b, 1e-5);
        
        lam = (gs_res.a + gs_res.b) / 2.0;
        L = 1.0 + 8.0 * lam;
        double mu = 1.0;
        double k_cond = L / mu;
        double h = 1.0 / L;
        
        double beta_hb = std::pow((std::sqrt(L) - std::sqrt(mu)) / (std::sqrt(L) + std::sqrt(mu)), 2);
        double beta_nag = (std::sqrt(L) - std::sqrt(mu)) / (std::sqrt(L) + std::sqrt(mu));
        
        std::printf("lambda* = %.6f\n", lam);
        std::printf("L = %.6f, kappa = %.6f, h = %.6f\n", L, k_cond, h);
        std::printf("Оптимальне beta (кулька) = %.6f\n", beta_hb);
        std::printf("Оптимальне beta (Нестеров) = %.6f\n\n", beta_nag);

        std::printf("Отримання еталонного розв'язку (довгий спуск 3000 ітерацій)...\n");
        GradResult ref_res = gradient_method(result.noisy, lam, L, "hb", beta_hb, h, 3000, 1e-12);
        Image I_star = ref_res.I;
        double E_star = ref_res.E;
        std::printf("Точне E* = %.8f (Норма градієнта: %.2e)\n\n", E_star, ref_res.grad);

        struct RunParams { std::string name, kind; double beta; };
        std::vector<RunParams> runs = {
            {"Спуск, h=1/L",             "gd",     0.0},
            {"Спуск, Армійо",            "armijo", 0.0},
            {"Важка кулька, beta_opt",   "hb",     beta_hb},
            {"Важка кулька, beta=0.9",   "hb",     0.9},
            {"Нестеров, опуклий",        "nag",    0.0},
            {"Нестеров, сильно опуклий", "nag",    beta_nag}
        };

        std::printf("%-26s | Ітер. | Nf  | Ngrad | E(I_N)\n", "Метод");
        std::printf("-----------------------------------------------------------------\n");
        
        for (const auto& run : runs) {
            GradResult res = gradient_method(result.noisy, lam, L, run.kind, run.beta, -1.0, 200, 1e-6, &I_star);
            std::printf("%-26s | %5d | %3ld | %5ld | %.8f %s\n", 
                run.name.c_str(), res.iterations, res.n_f, res.n_grad, res.E, 
                (res.converged ? "" : " (не збігся)"));
                
            std::string filename = "lab4_" + run.kind + "_" + std::to_string((int)(run.beta*100)) + ".txt";
            std::ofstream out(filename);
            for (const auto& r : res.rows) {
                out << r.k << " " << std::fixed << std::setprecision(10) << r.E << " " << r.dist << "\n";
            }
            out.close();
        }

        Image restored = I_star;
        for (size_t i = 0; i < restored.px.size(); ++i) {
            restored.px[i] = std::min(std::max(restored.px[i], 0.0), 1.0);
        }
        write_pgm(restored, "lab4_restored.pgm");

        std::printf("\n Відновлення зображення \n");
        std::printf("MSE  = %.8f\n", mse(restored, clean));
        std::printf("PSNR = %.4f дБ\n", psnr(restored, clean));

        std::printf("\n=== МАТЕМАТИЧНІ ПЕРЕВІРКИ ===\n");
        Image V(clean.h, clean.w);
        for (size_t i = 0; i < V.px.size(); ++i) V.px[i] = (rand() % 1000) / 1000.0;
        double t = 1e-5;
        
        Image I_t = clean;
        for (size_t i = 0; i < I_t.px.size(); ++i) I_t.px[i] += t * V.px[i];
        
        double E0 = energy(clean, result.noisy, lam);
        double E1 = energy(I_t, result.noisy, lam);
        double num_grad = (E1 - E0) / t;
        
        Image anal_grad = energy_grad(clean, result.noisy, lam);
        double dot_prod = 0.0;
        for (size_t i = 0; i < V.px.size(); ++i) dot_prod += anal_grad.px[i] * V.px[i];
        
        std::printf("Скінченна різниця dE/dt : %.6f\n", num_grad);
        std::printf("Аналітичний скаляр. добуток: %.6f\n", dot_prod);
        std::printf("Похибка градієнта       : %.2e\n", std::abs(num_grad - dot_prod));

        // Додана перевірка спряженості <grad(u), v> == <-u, div(v)>
        Field grad_u = grad_xy(clean);
        Field field_v = grad_xy(raw_img); 
        Image div_v = divergence(field_v.x, field_v.y);
        double left_conj = 0.0, right_conj = 0.0;
        for (size_t i = 0; i < clean.px.size(); ++i) {
            left_conj += grad_u.x.px[i] * field_v.x.px[i] + grad_u.y.px[i] * field_v.y.px[i];
            right_conj += clean.px[i] * (-div_v.px[i]);
        }
        std::printf("Спряженість: <grad(u), v> = %f, <-u, div(v)> = %f\n\n", left_conj, right_conj);

        std::printf("=== ДОСЛІД З ВЕЛИКИМ К (lambda = 50) ===\n");
        double lam_hard = 50.0;
        double L_hard = 1.0 + 8.0 * lam_hard;
        double k_cond_hard = L_hard / mu;
        double h_hard = 1.0 / L_hard;
        double beta_nag_hard = (std::sqrt(L_hard) - std::sqrt(mu)) / (std::sqrt(L_hard) + std::sqrt(mu));
        
        std::printf("kappa = %.2f\n", k_cond_hard);
        
        auto res_gd_hard = gradient_method(result.noisy, lam_hard, L_hard, "gd", 0.0, h_hard, 500, 1e-6);
        auto res_nag_hard = gradient_method(result.noisy, lam_hard, L_hard, "nag", beta_nag_hard, h_hard, 500, 1e-6);
        
        std::printf("Спуск (h=1/L) : %d ітер, E = %.6f\n", res_gd_hard.iterations, res_gd_hard.E);
        std::printf("Нестеров      : %d ітер, E = %.6f\n", res_nag_hard.iterations, res_nag_hard.E);
        
        std::ofstream out_gd_hard("lab4_gd_hard.txt");
        for (const auto& r : res_gd_hard.rows) out_gd_hard << r.k << " " << r.E << "\n";
        
        std::ofstream out_nag_hard("lab4_nag_hard.txt");
        for (const auto& r : res_nag_hard.rows) out_nag_hard << r.k << " " << r.E << "\n";
    }

    // лаб 5
    {
        std::printf("\n=== ЧАСТИНА 5: УМОВНА ОПТИМІЗАЦІЯ ===\n");
        Violation v_raw = get_violation(raw_img);
        
        std::printf(" Підхід A: Проєкція градієнта \n");
        PGResult res_pg = projected_gradient(raw_img, lam, L, false); 
        std::printf("Ітерацій: %d, Nf: %ld, Разом звернень: %ld\n", res_pg.iterations, res_pg.n_f, res_pg.n_f + res_pg.n_grad);
        std::printf("E* = %.8f, Макс. порушення: %.8f\n\n", res_pg.E, res_pg.rows.back().viol);

        std::printf(" Підхід A: FISTA (з прискоренням)\n");
        PGResult res_fista = projected_gradient(raw_img, lam, L, true); 
        std::printf("Ітерацій: %d, Nf: %ld, Разом звернень: %ld\n", res_fista.iterations, res_fista.n_f, res_fista.n_f + res_fista.n_grad);
        std::printf("E* = %.8f, Макс. порушення: %.8f\n\n", res_fista.E, res_fista.rows.back().viol);

        std::printf("Підхід B: Метод зовнішніх штрафів \n");
        PenaltyResult res_pen = penalty_method(raw_img, lam, L);
        std::printf("Зовнішніх ітерацій: %d, Nf: %ld, Разом звернень: %ld\n", (int)res_pen.rows.size(), res_pen.n_f, res_pen.n_f + res_pen.n_grad);
        std::printf("E* = %.8f, Порушення: %.2e\n\n", res_pen.E, res_pen.rows.back().violation);

        GradResult res_unc = gradient_method(raw_img, lam, L, "gd", 0.0, -1.0, 200, 1e-6);
        Violation v_unc = get_violation(res_unc.I);

        std::printf("=== ПІДСУМОК ЛАБ 5: ПІКСЕЛІ ПОЗА МЕЖАМИ [0, 1] ===\n");
        std::printf("Сире спостереження: %d px\n", v_raw.count);
        std::printf("Безумовний розв'язок: %d px\n", v_unc.count);
        std::printf("Проєкція градієнта: %d px\n", res_pg.rows.back().active);
        std::printf("Метод штрафів: %d px\n\n", res_pen.rows.back().outside);

        std::printf(" РЕЗУЛЬТАТ ВІДНОВЛЕННЯ (ПРОЄКЦІЯ) \n");
        std::printf("PSNR = %.2f дБ\n", psnr(res_pg.I, clean));
        std::printf("MSE  = %.3e\n", mse(res_pg.I, clean));

        write_pgm(res_pg.I, "lab5_restored.pgm");

        std::ofstream hist_raw("lab5_hist_raw.txt");
        for (double p : raw_img.px) hist_raw << p << "\n";
        hist_raw.close();

        std::ofstream hist_unc("lab5_hist_unc.txt");
        for (double p : res_unc.I.px) hist_unc << p << "\n";
        hist_unc.close();

        std::ofstream out_pg("lab5_pg.txt");
        for (const auto& r : res_pg.rows) out_pg << r.k << " " << std::fixed << std::setprecision(10) << r.E << "\n";
        out_pg.close();

        std::ofstream out_fista("lab5_fista.txt");
        for (const auto& r : res_fista.rows) out_fista << r.k << " " << std::fixed << std::setprecision(10) << r.E << "\n";
        out_fista.close();
        
        std::ofstream out_pen("lab5_pen.txt");
        int pen_k = 0;
        for (const auto& r : res_pen.rows) out_pen << pen_k++ << " " << std::fixed << std::setprecision(10) << r.E << "\n";
        out_pen.close();
    }
    return 0;
}
