/*
    실제 채널 SNR(TRUE_SNR_DB, 기본 3dB)은 고정한 채, 디코더가 LLR 계산에 쓰는
    sigma^2 추정용 SNR("assumed SNR")을 -3dB~+3dB 오차 범위로 스윕하며
    BER/FER이 어떻게 변하는지 관찰합니다.

    - 채널에 실제로 들어가는 잡음(sigma2_true)은 TRUE_SNR_DB로 고정.
    - Frozen Mask도 TRUE_SNR_DB(정확한 채널 조건)로 고정 설계 — 마스크 자체는 옳다고 가정.
    - 오직 LLR 계산에 쓰이는 sigma2_est만, "디코더가 착각한 SNR(assumed SNR
      = TRUE_SNR_DB + estimation_error_db)"로 계산해서 오차를 흉내낸다.
    - x축 = estimation_error_db (assumed - true), y축 = BER, FER (로그 스케일).
      error_db = 0인 지점이 곧 "정확하게 추정했을 때"의 결과.

    (예: 교수님 예시처럼 실제 3dB인데 2dB로 착각하면 estimation_error_db = 2 - 3 = -1dB)

    gcc -O2 -Wall -o polar_code_sigma_error_sweep polar_code_sigma_error_sweep.c -lm
    ./polar_code_sigma_error_sweep
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

// 실제 채널의 SNR (고정). 교수님 예시에서는 3dB.
#define TRUE_SNR_DB 3.0

// "디코더가 착각한 정도"를 이 범위로 스윕합니다 (단위: dB, assumed - true).
// 예를 들어 ERROR_DB = -1.0이면, 디코더는 실제보다 1dB 낮은 SNR로 sigma^2를 계산합니다.
#define ERROR_START -3.0
#define ERROR_END    3.0
#define ERROR_STEP   0.25

#define TARGET_ERRORS 1000

#define RESULT_DIR "result"
#define RESULT_FILE RESULT_DIR "/sigma_estimation_error_sweep.txt"

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
// 난수 생성기 (Xorshift64)
// =================================================================
static uint64_t rng_state;

static void seed_xorshift64(void) {
    rng_state = (uint64_t)time(NULL) ^ 0x5DEECE66DULL;
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

// Box-Muller 가우시안 난수 생성 (AWGN 채널용)
static double gauss_box_muller(void) {
    double u1, u2;
    do {
        u1 = rand_double();
    } while (u1 <= 1e-12);
    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// =================================================================
// Polar Code 구성 (Gaussian Approximation 기반 Frozen Mask)
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

// SNR(dB) -> sigma^2 (잡음 분산) 변환. 이 함수 하나로 "실제 SNR"과
// "디코더가 착각한(assumed) SNR"을 모두 sigma^2로 바꿔서 씁니다.
static double compute_sigma2_from_ebno_db(double ebno_db) {
    double ebno_linear = pow(10.0, ebno_db / 10.0);
    return 1.0 / (2.0 * R * ebno_linear);
}

// sigma2를 기준으로 각 서브채널의 신뢰도(mu)를 계산하고,
// 신뢰도가 높은 상위 K개 인덱스를 정보 비트 위치(info_mask=1)로 지정
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

    // bit-reversal을 통해 natural order로 변환
    double natural_mu[N];
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < n; j++) {
            if (i & (1 << j)) rev |= (1 << (n - 1 - j));
        }
        natural_mu[rev] = mu[i];
    }

    // 신뢰도 내림차순 정렬 (선택 정렬)
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

// Natural order 재귀형 Polar 인코더
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

// LLR 결합 함수 (f-function) — 이 실험에서는 일부러 Min-Sum이 아니라
// Exact SPA(log(1+exp(...)) 포함 정확한 로그 결합)를 사용합니다.
//
// 왜? Min-Sum(sign*min(|L1|,|L2|))과 g-function은 둘 다 "1차 동차(linear)"
// 함수라서, 모든 LLR에 똑같은 양수 배율(2/sigma2_est)이 곱해져도 재귀 전체에서
// 부호가 절대 바뀌지 않습니다. SC의 최종 결정은 LLR "부호"만 보므로, sigma2_est를
// 얼마로 착각하든(0보다 크기만 하면) 최종 비트 결정이 완전히 똑같아져
// BER/FER이 이론적으로 전혀 변하지 않는 평평한 직선이 되어버립니다.
// (직접 시뮬레이션으로 확인한 결과이기도 합니다.)
//
// 반면 Exact SPA의 log(1+exp(-(a1+a2))) - log(1+exp(-|a1-a2|)) 보정항은
// LLR의 "절대적인 크기"에 따라 비선형적으로 값이 달라지므로, sigma2_est
// 오차가 실제로 디코딩 결과에 영향을 줍니다 — 그래서 이 실험에는 이 방식을 씁니다.
static void polar_sc_decode_recursive(const double *llr, const int *info_mask,
                                      int *u_hat, int *u_coded,
                                      int n_len, int *mask_offset) {
    if (n_len == 1) {
        if (info_mask[*mask_offset]) {
            u_hat[0] = (llr[0] < 0.0) ? 1 : 0;
        } else {
            u_hat[0] = 0; // frozen bit
        }
        u_coded[0] = u_hat[0];
        (*mask_offset)++;
        return;
    }

    int n2 = n_len / 2;
    double llr_left[N / 2], llr_right[N / 2];
    int u_hat_left[N / 2], u_hat_right[N / 2];
    int u_coded_left[N / 2], u_coded_right[N / 2];

    // f-function (Exact SPA)
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

    // g-function (결정 피드백)
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

// frozen mask에 따라 정보 비트를 u 벡터에 배치 (frozen 위치는 0)
static void build_u_from_mask(const int *mask, const int *info_bits, int *u) {
    int info_ptr = 0;
    for (int i = 0; i < N; i++) {
        if (mask[i]) u[i] = info_bits[info_ptr++];
        else u[i] = 0;
    }
}

// BPSK 변조 + AWGN 채널 (실제 채널이므로 항상 sigma2_true로 잡음을 생성)
static void awgn_channel(const int *x, double sigma2_true, double *rx) {
    double sigma = sqrt(sigma2_true);
    for (int i = 0; i < N; i++) {
        double tx = 1.0 - 2.0 * x[i]; // 0 -> +1, 1 -> -1
        rx[i] = tx + gauss_box_muller() * sigma;
    }
}

// LLR 계산 — 여기서 쓰는 sigma2_est가 "디코더가 착각한 잡음 세기"입니다.
// sigma2_est == sigma2_true면 정확한 추정, 다르면 그만큼 오차가 있는 것.
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

typedef struct {
    double ber;
    double fer;
} BerFerResult;

// 한 조건(sigma2_true=실제 채널, sigma2_est=디코더가 착각한 값)에 대해
// 오류가 TARGET_ERRORS개 모일 때까지 몬테카를로 반복하는 코어 함수.
static BerFerResult run_monte_carlo(const int *info_mask, double sigma2_true, double sigma2_est) {
    long total_errors = 0, total_bits = 0, total_frames = 0, total_frame_errors = 0;

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

    BerFerResult result;
    result.ber = (double)total_errors / (double)total_bits;
    result.fer = (double)total_frame_errors / (double)total_frames;
    return result;
}

// =================================================================
// 시뮬레이션 — sigma^2 추정 오차(estimation_error_db) 스윕
// =================================================================
static void run_error_sweep(FILE *fp) {
    // 실제 채널은 TRUE_SNR_DB로 고정. Frozen Mask도 "실제 채널을 정확히 안다"는
    // 가정으로 TRUE_SNR_DB 기준으로 한 번만 설계하고, 오차 스윕 내내 그대로 사용한다.
    double sigma2_true = compute_sigma2_from_ebno_db(TRUE_SNR_DB);
    int info_mask[N];
    construct_frozen_mask(sigma2_true, info_mask);

    printf("\nTrue channel SNR = %.2f dB (fixed)\n", TRUE_SNR_DB);
    printf("Estimation Error (dB)\tAssumed SNR (dB)\tBER\t\tFER\n");

    for (double error_db = ERROR_START; error_db <= ERROR_END + 1e-9; error_db += ERROR_STEP) {
        double assumed_snr_db = TRUE_SNR_DB + error_db;
        double sigma2_est = compute_sigma2_from_ebno_db(assumed_snr_db);

        BerFerResult result = run_monte_carlo(info_mask, sigma2_true, sigma2_est);

        printf("%+.2f\t\t\t%.2f\t\t\t%.6e\t%.6e\n",
               error_db, assumed_snr_db, result.ber, result.fer);
        fprintf(fp, "%.2f %.2f %.6e %.6e\n", error_db, assumed_snr_db, result.ber, result.fer);
    }
}

// gnuplot으로 "SNR 추정 오차 vs BER/FER" 그래프를 그린다.
// x=0(오차 없음, 정확한 추정) 지점에 기준선을 표시한다.
static void run_gnuplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Warning: failed to launch gnuplot. Is it installed and in PATH?\n");
        return;
    }

    fprintf(gp, "set terminal windows size 900, 650\n");
    fprintf(gp, "set title 'Polar Code (N=1024, K=512) SC Decoder: sigma^2 Estimation Error "
                "(True SNR = %.1f dB fixed)'\n", TRUE_SNR_DB);
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'SNR Estimation Error (dB) = Assumed SNR - True SNR'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp, "set arrow from 0, graph 0 to 0, graph 1 nohead lc rgb 'gray' dt 2\n");
    fprintf(gp, "set label 'error = 0 (정확한 추정)' at 0, graph 0.95 offset 1,0 tc rgb 'gray'\n");
    fprintf(gp,
            "plot '" RESULT_FILE "' using 1:3 with linespoints lw 2 pt 7 lc rgb 'red' title 'BER', "
            "'" RESULT_FILE "' using 1:4 with linespoints lw 2 pt 7 lc rgb 'blue' title 'FER'\n");

    PCLOSE(gp);
}

int main(void) {
    seed_xorshift64();
    ensure_result_dir();

    FILE *fp = fopen(RESULT_FILE, "w");
    if (fp == NULL) {
        printf("Failed to open output file: %s\n", RESULT_FILE);
        return 1;
    }

    printf("Running Polar Code (N=%d, K=%d) SC simulation: sigma^2 estimation error sweep...\n", N, K);
    run_error_sweep(fp);
    fclose(fp);

    printf("\nSimulation completed! Results saved to %s\n", RESULT_FILE);

    run_gnuplot();
    printf("Gnuplot window launched (if gnuplot is installed).\n");

    return 0;
}