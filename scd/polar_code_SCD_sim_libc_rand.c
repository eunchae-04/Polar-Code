/*
    Xorshift64 대신 표준 라이브러리 rand()/srand()를 사용하는 대조군(control) 버전입니다.
    난수 생성기 부분(seed_xorshift64, xorshift64, rand_double)만 바뀌고,
    나머지 알고리즘(Frozen Mask, 인코더, SC 디코더, 시뮬레이션 루프)은 원본과 완전히 동일합니다.
    이 결과가 polar_code_SCD_sim.c(Xorshift64 버전) 결과와 비슷하게 나오면,
    "커스텀 RNG를 써도 결과가 왜곡되지 않는다"는 걸 검증할 수 있습니다.

    gcc -O2 -Wall -o polar_code_SCD_sim_libc_rand polar_code_SCD_sim_libc_rand.c -lm
    ./polar_code_SCD_sim_libc_rand
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>

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

#define RESULT_DIR "result"
#define RESULT_FILE RESULT_DIR "/simulation_results_libc_rand.txt"

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(path) mkdir(path, 0755)
#endif

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

static void ensure_dir_recursive(const char *path) {
    char buffer[256];
    size_t length = strlen(path);

    if (length >= sizeof(buffer)) {
        printf("Warning: output path too long: %s\n", path);
        return;
    }

    strcpy(buffer, path);

    for (char *cursor = buffer + 1; *cursor; ++cursor) {
        if (*cursor == '/' || *cursor == '\\') {
            char saved = *cursor;
            *cursor = '\0';
            if (MKDIR(buffer) != 0 && errno != EEXIST) {
                printf("Warning: failed to create '%s' directory (errno=%d). Output file may fail to save.\n",
                       buffer, errno);
                *cursor = saved;
                return;
            }
            *cursor = saved;
        }
    }

    if (MKDIR(buffer) != 0 && errno != EEXIST) {
        printf("Warning: failed to create '%s' directory (errno=%d). Output file may fail to save.\n",
               buffer, errno);
    }
}

static void ensure_result_dir(void) {
    ensure_dir_recursive(RESULT_DIR);
}

// =================================================================
// 난수 생성기 — 여기만 표준 라이브러리 rand()/srand()로 교체
// (원본의 seed_xorshift64/xorshift64 자리를 대신함)
// =================================================================
static void seed_libc_rand(void) {
    srand((unsigned int)time(NULL));
}

// xorshift64()는 64비트 정수를 돌려줬지만, rand()는 통상 0~32767(최소 보장 범위)
// 밖에 못 돌려주므로, 비트를 여러 번 뽑아 이어붙여 uint64_t 하나를 만들어준다.
// (rand_double, generate_info_bits 등 기존 코드가 요구하는 인터페이스를 그대로 맞추기 위함)
static uint64_t libc_rand64(void) {
    uint64_t r = 0;
    for (int i = 0; i < 5; i++) {          // RAND_MAX가 최소 15비트이므로 5번이면 64비트 이상 확보
        r = (r << 15) ^ (uint64_t)(rand() & 0x7FFF);
    }
    return r;
}

static double rand_double(void) {
    return (double)libc_rand64() / (double)UINT64_MAX;
}

// Box-Muller 가우시안 난수 생성 (AWGN 채널용) — 원본과 완전히 동일
static double gauss_box_muller(void) {
    double u1, u2;
    do {
        u1 = rand_double();
    } while (u1 <= 1e-12);
    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// =================================================================
// Polar Code 구성 (Gaussian Approximation 기반 Frozen Mask) — 원본과 동일
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

    int n = 10;
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
// Encode / Decode — 원본과 동일
// =================================================================
static void polar_encode_recursive(const int *u, int *x, int n_len) {
    if (n_len == 1) {
        x[0] = u[0];
        return;
    }

    int n2 = n_len / 2;
    int v[N / 2];

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
    double llr_left[N / 2], llr_right[N / 2];
    int u_hat_left[N / 2], u_hat_right[N / 2];
    int u_coded_left[N / 2], u_coded_right[N / 2];

    for (int i = 0; i < n2; i++) {
        double l1 = llr[i];
        double l2 = llr[i + n2];
        double sign = ((l1 < 0.0) ^ (l2 < 0.0)) ? -1.0 : 1.0;
        double a1 = fabs(l1);
        double a2 = fabs(l2);
        double min_val = (a1 < a2) ? a1 : a2;
        llr_left[i] = sign * min_val;
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
// 공통 유틸리티 — generate_info_bits만 libc_rand64() 사용, 나머지는 동일
// =================================================================
static void generate_info_bits(int *info_bits) {
    for (int i = 0; i < K; i++) {
        info_bits[i] = (int)(libc_rand64() & 1ULL);
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

// =================================================================
// 시뮬레이션 (Eb/No 스윕 + 몬테카를로) — 원본과 동일
// =================================================================
static void run_simulation(FILE *fp) {
    printf("\nEb/No (dB)\tBER\t\tFER\n");

    for (double snr_db = EBNO_START; snr_db <= EBNO_END + 1e-9; snr_db += EBNO_STEP) {
        double sigma2_true = compute_sigma2_from_ebno_db(snr_db);
        double sigma2_est = sigma2_true;

        int info_mask[N];
        construct_frozen_mask(sigma2_true, info_mask);

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
}

static void run_gnuplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Warning: failed to launch gnuplot. Is it installed and in PATH?\n");
        return;
    }

    fprintf(gp, "set terminal windows size 800, 600\n");
    fprintf(gp, "set title 'Polar Code SC Decoder Performance (libc rand() control group)'\n");
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp,
            "plot '" RESULT_FILE "' using 1:2 with linespoints lw 2 lc rgb 'red' title 'BER', "
            "'" RESULT_FILE "' using 1:3 with linespoints lw 2 lc rgb 'blue' title 'FER'\n");

    PCLOSE(gp);
}

int main(void) {
    seed_libc_rand();
    ensure_result_dir();

    FILE *fp = fopen(RESULT_FILE, "w");
    if (fp == NULL) {
        printf("Failed to open output file: %s\n", RESULT_FILE);
        return 1;
    }

    printf("Running Polar Code (N=%d, K=%d) SC simulation (libc rand() control group)...\n", N, K);
    run_simulation(fp);
    fclose(fp);

    printf("\nSimulation completed! Results saved to %s\n", RESULT_FILE);

    run_gnuplot();
    printf("Gnuplot window launched (if gnuplot is installed).\n");

    return 0;
}