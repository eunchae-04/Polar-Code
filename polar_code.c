#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 1024
#define K 512
#define R ((double)K / (double)N)

#define EBNO_START 0.0
#define EBNO_END   3.0
#define EBNO_STEP  0.25
#define TARGET_ERRORS 10000

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

typedef struct {
    const char *output_file;
    const char *title;
    double llr_scale;             // sigma2_est = sigma2_true * llr_scale
    bool use_fixed_mask;          // true면 fixed_mask_snr_db로 생성한 마스크를 고정 사용
    double fixed_mask_snr_db;     // fixed mask 생성용 SNR
} SimulationConfig;

// =================================================================
// 난수 생성기
// =================================================================
static uint64_t rng_state;

static void seed_xorshift64(void) {
    rng_state = (uint64_t)time(NULL) ^ 0x5DEECE66DL;
    if (rng_state == 0) rng_state = 1;
}

static uint64_t xorshift64(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static double rand_double(void) {
    return (double)xorshift64() / (double)UINT64_MAX;
}

static double gauss_box_muller(void) {
    double u1, u2;
    do {
        u1 = rand_double();
    } while (u1 <= 1e-12);
    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// =================================================================
// Polar Code 구성
// =================================================================
static double phi(double x) {
    if (x <= 0.0) return 1.0;
    if (x <= 10.0) return exp(-0.4527 * pow(x, 0.86) + 0.0218);
    return sqrt(M_PI / x) * (1.0 - 10.0 / (7.0 * x)) * exp(-x / 4.0);
}

static double phi_inv(double y) {
    if (y >= 1.0) return 0.0;
    if (y < 1e-200) y = 1e-200;

    double low = 1e-6;
    double high = 1000.0;
    double mid = low;

    for (int iter = 0; iter < 35; iter++) {
        mid = (low + high) / 2.0;
        if (phi(mid) > y) low = mid;
        else high = mid;
    }
    return mid;
}

static double compute_sigma2_from_ebno_db(double ebno_db) {
    double ebno_linear = pow(10.0, ebno_db / 10.0);
    return 1.0 / (2.0 * R * ebno_linear);
}

static void construct_frozen_mask(double sigma2, int *info_mask) {
    double mu[N];
    mu[0] = 2.0 / sigma2;

    int n = 10; // N = 2^10 = 1024
    for (int stage = 1; stage <= n; stage++) {
        int block_size = 1 << (stage - 1);
        for (int i = 0; i < block_size; i++) {
            double T = mu[i];
            mu[i] = phi_inv(1.0 - pow(1.0 - phi(T), 2.0));
            mu[i + block_size] = 2.0 * T;
        }
    }

    double natural_mu[N];
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < n; j++) {
            if (i & (1 << j)) rev |= (1 << (n - 1 - j));
        }
        natural_mu[rev] = mu[i];
    }

    int idx[N];
    for (int i = 0; i < N; i++) idx[i] = i;

    for (int i = 0; i < N - 1; i++) {
        int max_j = i;
        for (int j = i + 1; j < N; j++) {
            if (natural_mu[idx[j]] > natural_mu[idx[max_j]]) max_j = j;
        }
        int tmp = idx[i];
        idx[i] = idx[max_j];
        idx[max_j] = tmp;
    }

    for (int i = 0; i < N; i++) info_mask[i] = 0;
    for (int i = 0; i < K; i++) info_mask[idx[i]] = 1;
}

// =================================================================
// Encode / Decode
// =================================================================
static void polar_encode_recursive(const int *u, int *x, int n_len) {
    if (n_len == 1) {
        x[0] = u[0];
        return;
    }

    int n2 = n_len / 2;
    int v[n2];

    for (int i = 0; i < n2; i++) {
        v[i] = u[i] ^ u[i + n2];
    }

    polar_encode_recursive(v, x, n2);
    polar_encode_recursive(u + n2, x + n2, n2);
}

static void polar_sc_decode_recursive(const double *llr, const int *info_mask,
                                      int *u_hat, int *u_coded,
                                      int n_len, int *mask_offset) {
    if (n_len == 1) {
        if (info_mask[*mask_offset]) {
            u_hat[0] = (llr[0] < 0.0) ? 1 : 0;
        } else {
            u_hat[0] = 0;
        }
        u_coded[0] = u_hat[0];
        (*mask_offset)++;
        return;
    }

    int n2 = n_len / 2;
    double llr_left[n2], llr_right[n2];
    int u_hat_left[n2], u_hat_right[n2];
    int u_coded_left[n2], u_coded_right[n2];

    for (int i = 0; i < n2; i++) {
        double l1 = llr[i];
        double l2 = llr[i + n2];
        double sign = ((l1 < 0.0) ^ (l2 < 0.0)) ? -1.0 : 1.0;
        double a1 = fabs(l1);
        double a2 = fabs(l2);
        double min_val = (a1 < a2) ? a1 : a2;

        llr_left[i] = sign * min_val
                    + log(1.0 + exp(-(a1 + a2)))
                    - log(1.0 + exp(-fabs(a1 - a2)));
    }

    polar_sc_decode_recursive(llr_left, info_mask, u_hat_left, u_coded_left, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        llr_right[i] = llr[i + n2] + (1 - 2 * u_coded_left[i]) * llr[i];
    }

    polar_sc_decode_recursive(llr_right, info_mask, u_hat_right, u_coded_right, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        u_hat[i] = u_hat_left[i];
        u_hat[i + n2] = u_hat_right[i];
        u_coded[i] = u_coded_left[i] ^ u_coded_right[i];
        u_coded[i + n2] = u_coded_right[i];
    }
}

// =================================================================
// 공통 유틸리티
// =================================================================
static void generate_info_bits(int *info_bits) {
    for (int i = 0; i < K; i++) {
        info_bits[i] = (int)(xorshift64() & 1ULL);
    }
}

static void build_u_from_mask(const int *mask, const int *info_bits, int *u) {
    int info_ptr = 0;
    for (int i = 0; i < N; i++) {
        if (mask[i]) u[i] = info_bits[info_ptr++];
        else u[i] = 0;
    }
}

static void awgn_channel(const int *x, double sigma2_true, double *rx) {
    double sigma = sqrt(sigma2_true);
    for (int i = 0; i < N; i++) {
        double tx = 1.0 - 2.0 * x[i];
        rx[i] = tx + gauss_box_muller() * sigma;
    }
}

static void compute_llr(const double *rx, double sigma2_est, double *llr) {
    for (int i = 0; i < N; i++) {
        llr[i] = (2.0 / sigma2_est) * rx[i];
    }
}

static long count_bit_errors(const int *info_mask, const int *u_hat, const int *info_bits) {
    int info_ptr = 0;
    long bit_errors = 0;

    for (int i = 0; i < N; i++) {
        if (info_mask[i]) {
            if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
        }
    }
    return bit_errors;
}

static void run_simulation(const SimulationConfig *cfg, const int *fixed_mask) {
    FILE *fp = fopen(cfg->output_file, "w");
    if (fp == NULL) {
        printf("Failed to open output file: %s\n", cfg->output_file);
        return;
    }

    printf("\n%s\n", cfg->title);
    printf("Eb/No (dB)\tBER\t\tFER\n");

    for (double snr_db = EBNO_START; snr_db <= EBNO_END + 1e-9; snr_db += EBNO_STEP) {
        double sigma2_true = compute_sigma2_from_ebno_db(snr_db);
        double sigma2_est = sigma2_true * cfg->llr_scale;

        int info_mask[N];
        if (cfg->use_fixed_mask) {
            for (int i = 0; i < N; i++) info_mask[i] = fixed_mask[i];
        } else {
            construct_frozen_mask(sigma2_true, info_mask);
        }

        long total_errors = 0;
        long total_bits = 0;
        long total_frames = 0;
        long total_frame_errors = 0;

        while (total_errors < TARGET_ERRORS) {
            int info_bits[K];
            int u[N];
            int x[N];
            double rx[N];
            double llr[N];
            int u_hat[N];
            int u_coded[N];
            int mask_offset = 0;

            generate_info_bits(info_bits);
            build_u_from_mask(info_mask, info_bits, u);
            polar_encode_recursive(u, x, N);
            awgn_channel(x, sigma2_true, rx);
            compute_llr(rx, sigma2_est, llr);
            polar_sc_decode_recursive(llr, info_mask, u_hat, u_coded, N, &mask_offset);

            long bit_errors = count_bit_errors(info_mask, u_hat, info_bits);

            if (bit_errors > 0) total_frame_errors++;
            total_errors += bit_errors;
            total_bits += K;
            total_frames++;
        }

        double final_ber = (double)total_errors / (double)total_bits;
        double final_fer = (double)total_frame_errors / (double)total_frames;

        printf("%.2f dB\t\t%.6e\t%.6e\n", snr_db, final_ber, final_fer);
        fprintf(fp, "%.2f %.6e %.6e\n", snr_db, final_ber, final_fer);
    }

    fclose(fp);
}

static void run_gnuplot_multiplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) return;

    fprintf(gp, "set terminal windows size 1300, 450\n");
    fprintf(gp, "set multiplot layout 1,3 title 'Polar Code Scenario Analyses (N=1024, K=512)' font ',13'\n");
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "unset yrange\n");
    fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");

    fprintf(gp, "set title '1. Baseline (Perfect)'\n");
    fprintf(gp,
            "plot 'simulation_baseline.txt' using 1:2 with linespoints lw 2 lc rgb 'purple' title 'BER', "
            "'simulation_baseline.txt' using 1:3 with linespoints lw 2 lc rgb 'cyan' title 'FER'\n");

    fprintf(gp, "set title '2. Scenario A (LLR Error, a=0.5)'\n");
    fprintf(gp,
            "plot 'simulation_scenario_a.txt' using 1:2 with linespoints lw 2 lc rgb 'red' title 'BER', "
            "'simulation_scenario_a.txt' using 1:3 with linespoints lw 2 lc rgb 'orange' title 'FER'\n");

    fprintf(gp, "set title '3. Scenario B (Design SNR=0dB)'\n");
    fprintf(gp,
            "plot 'simulation_scenario_b.txt' using 1:2 with linespoints lw 2 lc rgb 'blue' title 'BER', "
            "'simulation_scenario_b.txt' using 1:3 with linespoints lw 2 lc rgb 'dark-green' title 'FER'\n");

    fprintf(gp, "unset multiplot\n");
    PCLOSE(gp);
}

int main(void) {
    seed_xorshift64();

    SimulationConfig baseline = {
        .output_file = "simulation_baseline.txt",
        .title = "[1/3] Running Baseline Scenario...",
        .llr_scale = 1.0,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_a = {
        .output_file = "simulation_scenario_a.txt",
        .title = "[2/3] Running Scenario A (LLR Mismatch, alpha = 0.5)...",
        .llr_scale = 0.5,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_b = {
        .output_file = "simulation_scenario_b.txt",
        .title = "[3/3] Running Scenario B (Design SNR Mismatch = 0.0 dB Fixed)...",
        .llr_scale = 1.0,
        .use_fixed_mask = true,
        .fixed_mask_snr_db = 0.0
    };

    int fixed_info_mask[N];
    double sigma2_design = compute_sigma2_from_ebno_db(scenario_b.fixed_mask_snr_db);
    construct_frozen_mask(sigma2_design, fixed_info_mask);

    run_simulation(&baseline, NULL);
    run_simulation(&scenario_a, NULL);
    run_simulation(&scenario_b, fixed_info_mask);

    run_gnuplot_multiplot();

    printf("\nAll simulations completed! Gnuplot multiplot window popped up.\n");
    return 0;
}