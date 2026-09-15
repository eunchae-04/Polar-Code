/*
    실제 채널 SNR(TRUE_SNR_DB, 기본 3dB)은 고정한 채, 디코더가 LLR 계산에 쓰는
    sigma^2 추정용 SNR("assumed SNR")을 -3dB~+3dB 오차 범위로 스윕하며
    BER/FER이 어떻게 변하는지 관찰합니다.

    이번 버전에서는 두 가지를 함께 봅니다.
      1) "정확하게 추정했을 때(error=0)"를 baseline으로 명확히 표시하고,
         그 대비 다른 오차 지점들이 몇 배 나빠지는지 숫자와 그래프로 비교.
      2) Min-Sum(선형, scale-invariant) 방식과 Exact SPA(비선형) 방식을
         같은 조건으로 동시에 시뮬레이션해서 나란히 비교.
         -> Min-Sum은 sigma^2를 얼마로 착각하든 최종 결정(u_hat)이 이론적으로
            바뀌지 않으므로 거의 평평한 직선이 나오고, Exact SPA는 오차에
            민감하게 반응해 뚜렷한 곡선이 나오는 것을 직접 확인할 수 있습니다.

    통계적 노이즈를 줄이기 위해 "공통 난수(Common Random Numbers)" 기법을 씁니다:
    모든 (오차값, 디코더) 조합이 정확히 같은 프레임 수 · 같은 시드로 시작해서
    완전히 동일한 정보비트·채널잡음 시퀀스를 겪습니다. 그래서 결과 차이는
    오직 sigma2_est(추정 오차) 하나 때문이라는 게 보장되고, 서로 다른 난수가
    섞여서 생기는 흔들림 없이 훨씬 적은 프레임 수로도 매끈한 곡선을 볼 수 있습니다.

    - 채널에 실제로 들어가는 잡음(sigma2_true)은 TRUE_SNR_DB로 고정.
    - Frozen Mask도 TRUE_SNR_DB(정확한 채널 조건)로 고정 설계 — 마스크 자체는 옳다고 가정.
    - 오직 LLR 계산에 쓰이는 sigma2_est만, "디코더가 착각한 SNR(assumed SNR
      = TRUE_SNR_DB + estimation_error_db)"로 계산해서 오차를 흉내낸다.

    (예: 실제 3dB인데 2dB로 착각하면 estimation_error_db = 2 - 3 = -1dB)

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

// 실제 채널의 SNR (고정). 예시: 3dB.
#define TRUE_SNR_DB 3.0

// "디코더가 착각한 정도"를 이 범위로 스윕합니다 (단위: dB, assumed - true).
#define ERROR_START -3.0
#define ERROR_END    3.0
#define ERROR_STEP   0.25

// (ERROR_END - ERROR_START) / ERROR_STEP + 1 = 25
// 범위/간격을 바꾸면 이 값도 함께 맞춰줘야 합니다.
#define NUM_ERROR_POINTS 25

// TARGET_ERRORS 기반 가변 반복 대신, "공통 난수(Common Random Numbers)" 기법을 씁니다:
// 모든 (오차값, 디코더) 조합에서 정확히 같은 프레임 수를, 정확히 같은 난수 시드로
// 시작해서 돌립니다. 그러면 모든 조건이 "완전히 동일한 정보비트·채널잡음 시퀀스"를
// 겪게 되어, 결과 차이가 오직 sigma2_est(추정 오차) 하나 때문이라는 게 보장됩니다.
// (서로 다른 무작위 잡음이 섞여서 생기는 노이즈가 사라지므로, TARGET_ERRORS 방식보다
//  훨씬 적은 프레임 수로도 매끈한 U자형 곡선을 볼 수 있습니다.)
#define NUM_FRAMES 430000

// 임의로 고른 고정 시드 (파이(π)의 16진수 앞자리). 시간 기반이 아니라 상수로 고정해서
// 모든 조건이 항상 같은 난수 시퀀스로 시작하도록 합니다.
#define COMMON_SEED 0x243F6A8885A308D3ULL

#define RESULT_DIR "result"
#define MINSUM_FILE   RESULT_DIR "/sigma_error_sweep_minsum.txt"
#define SPA_FILE      RESULT_DIR "/sigma_error_sweep_spa.txt"
#define SUMMARY_FILE  RESULT_DIR "/sigma_error_comparison_summary.txt"
#define BASELINE_MINSUM_FILE RESULT_DIR "/sigma_error_baseline_minsum.txt"
#define BASELINE_SPA_FILE    RESULT_DIR "/sigma_error_baseline_spa.txt"

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

// (참고: 이 실험은 공통 난수(Common Random Numbers) 기법을 쓰기 때문에
//  rng_state를 매번 COMMON_SEED로 고정 초기화합니다. 그래서 별도의
//  "실행할 때마다 다른 시드" 초기화 함수는 이 파일에서는 쓰지 않습니다.)

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

// SNR(dB) -> sigma^2 (잡음 분산) 변환.
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
// Encode
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

// =================================================================
// Decode — Min-Sum과 Exact SPA 두 가지 f-function을 모두 구현.
// g-function과 leaf(하드 디시전) 로직은 두 방식이 완전히 동일합니다.
// =================================================================

// Min-Sum: f-function = sign * min(|l1|, |l2|)  (선형, 1차 동차)
// -> LLR 전체에 같은 양수 배율(2/sigma2_est)이 곱해져도 부호가 안 바뀌므로,
//    sigma2_est를 얼마로 착각하든 최종 결정이 이론적으로 동일합니다.
static void polar_sc_decode_recursive_minsum(const double *llr, const int *info_mask,
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

    polar_sc_decode_recursive_minsum(llr_left, info_mask, u_hat_left, u_coded_left, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        llr_right[i] = llr[i + n2] + (1 - 2 * u_coded_left[i]) * llr[i];
    }

    polar_sc_decode_recursive_minsum(llr_right, info_mask, u_hat_right, u_coded_right, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        u_hat[i] = u_hat_left[i];
        u_hat[i + n2] = u_hat_right[i];
        u_coded[i] = u_coded_left[i] ^ u_coded_right[i];
        u_coded[i + n2] = u_coded_right[i];
    }
}

// Exact SPA: f-function = sign*min(|l1|,|l2|) + log(1+e^-|l1+l2|) - log(1+e^-|l1-l2|)
// (주의: 보정항은 반드시 부호가 있는 l1+l2, l1-l2의 절댓값을 써야 합니다.
//  |l1|+|l2|, ||l1|-|l2||로 계산하면 l1과 l2의 부호가 다를 때 두 항이 뒤바뀌어
//  틀린 값이 나옵니다 — 부호가 같을 때만 우연히 두 표현이 일치합니다.)
// -> 보정항이 절댓값 크기에 비선형적으로 반응하므로, sigma2_est 오차(=LLR 전체 배율 오차)가
//    실제로 결합 결과와 최종 결정에 영향을 줍니다.
static void polar_sc_decode_recursive_spa(const double *llr, const int *info_mask,
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
        llr_left[i] = sign * min_val
                    + log(1.0 + exp(-fabs(l1 + l2)))
                    - log(1.0 + exp(-fabs(l1 - l2)));
    }

    polar_sc_decode_recursive_spa(llr_left, info_mask, u_hat_left, u_coded_left, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        llr_right[i] = llr[i + n2] + (1 - 2 * u_coded_left[i]) * llr[i];
    }

    polar_sc_decode_recursive_spa(llr_right, info_mask, u_hat_right, u_coded_right, n2, mask_offset);

    for (int i = 0; i < n2; i++) {
        u_hat[i] = u_hat_left[i];
        u_hat[i + n2] = u_hat_right[i];
        u_coded[i] = u_coded_left[i] ^ u_coded_right[i];
        u_coded[i + n2] = u_coded_right[i];
    }
}

// 위 두 디코더를 같은 형태로 호출할 수 있도록 함수 포인터 타입을 정의.
typedef void (*DecodeFn)(const double *llr, const int *info_mask,
                          int *u_hat, int *u_coded, int n_len, int *mask_offset);

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

typedef struct {
    double ber;
    double fer;
} BerFerResult;

// 한 조건(sigma2_true=실제 채널, sigma2_est=디코더가 착각한 값, decode=사용할 디코더)에
// 대해 NUM_FRAMES개의 프레임을 처리한다. 매번 COMMON_SEED로 초기화하므로, 모든 조건이
// 완전히 같은 정보비트·채널잡음 시퀀스를 겪는다 (공통 난수 기법 — 노이즈 제거용).
static BerFerResult run_monte_carlo(const int *info_mask, double sigma2_true,
                                     double sigma2_est, DecodeFn decode) {
    rng_state = COMMON_SEED; // 모든 조건이 동일한 난수 시퀀스로 시작하도록 매번 리셋

    long total_errors = 0, total_bits = 0, total_frames = 0, total_frame_errors = 0;

    for (long frame = 0; frame < NUM_FRAMES; frame++) {
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
        decode(llr, info_mask, u_hat, u_coded, N, &mask_offset);

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
// Min-Sum과 Exact SPA를 같은 조건에서 함께 계산해서 비교한다.
// =================================================================
static void run_error_sweep(void) {
    double sigma2_true = compute_sigma2_from_ebno_db(TRUE_SNR_DB);
    int info_mask[N];
    construct_frozen_mask(sigma2_true, info_mask);

    double error_arr[NUM_ERROR_POINTS];
    double minsum_ber[NUM_ERROR_POINTS], minsum_fer[NUM_ERROR_POINTS];
    double spa_ber[NUM_ERROR_POINTS], spa_fer[NUM_ERROR_POINTS];
    int count = 0;

    printf("\nTrue channel SNR = %.2f dB (fixed), NUM_FRAMES = %d, total points = %d\n",
           TRUE_SNR_DB, NUM_FRAMES, NUM_ERROR_POINTS);
    printf("Error(dB)  MinSum_BER    MinSum_FER    SPA_BER       SPA_FER       "
           "| elapsed  ETA remaining\n");

    time_t t_start = time(NULL);

    for (double error_db = ERROR_START; error_db <= ERROR_END + 1e-9 && count < NUM_ERROR_POINTS; error_db += ERROR_STEP) {
        double assumed_snr_db = TRUE_SNR_DB + error_db;
        double sigma2_est = compute_sigma2_from_ebno_db(assumed_snr_db);

        BerFerResult r_minsum = run_monte_carlo(info_mask, sigma2_true, sigma2_est, polar_sc_decode_recursive_minsum);
        BerFerResult r_spa    = run_monte_carlo(info_mask, sigma2_true, sigma2_est, polar_sc_decode_recursive_spa);

        error_arr[count] = error_db;
        minsum_ber[count] = r_minsum.ber;
        minsum_fer[count] = r_minsum.fer;
        spa_ber[count] = r_spa.ber;
        spa_fer[count] = r_spa.fer;
        count++;

        // 진행 상황 + 남은 예상 시간(ETA) 계산 및 출력.
        // (지금까지 걸린 시간 / 완료한 지점 수) x 남은 지점 수 로 추정.
        double elapsed_sec = difftime(time(NULL), t_start);
        double sec_per_point = elapsed_sec / (double)count;
        double remaining_sec = sec_per_point * (double)(NUM_ERROR_POINTS - count);

        printf("%+6.2f     %.6e   %.6e   %.6e   %.6e   | %6.0fs  ~%.0fs left (%d/%d)\n",
               error_db, r_minsum.ber, r_minsum.fer, r_spa.ber, r_spa.fer,
               elapsed_sec, remaining_sec, count, NUM_ERROR_POINTS);
        fflush(stdout); // 콘솔에 즉시 출력되도록 강제 flush (버퍼링 때문에 안 보이는 것 방지)
    }

    // error=0(정확한 추정)을 baseline으로 찾는다. ERROR_STEP 격자에 0.00이 정확히
    // 포함되도록 범위를 잡았으므로, 부동소수점 오차만 감안해서 비교한다.
    int baseline_idx = -1;
    for (int i = 0; i < count; i++) {
        if (fabs(error_arr[i]) < ERROR_STEP / 2.0) { baseline_idx = i; break; }
    }

    // ---- 데이터 파일 저장 (gnuplot용) ----
    FILE *fp_minsum = fopen(MINSUM_FILE, "w");
    FILE *fp_spa = fopen(SPA_FILE, "w");
    for (int i = 0; i < count; i++) {
        if (fp_minsum) fprintf(fp_minsum, "%.2f %.6e %.6e\n", error_arr[i], minsum_ber[i], minsum_fer[i]);
        if (fp_spa) fprintf(fp_spa, "%.2f %.6e %.6e\n", error_arr[i], spa_ber[i], spa_fer[i]);
    }
    if (fp_minsum) fclose(fp_minsum);
    if (fp_spa) fclose(fp_spa);

    // baseline 지점만 담은 파일 (그래프에서 큰 별 마커로 강조하기 위함)
    if (baseline_idx >= 0) {
        FILE *fp_bm = fopen(BASELINE_MINSUM_FILE, "w");
        if (fp_bm) {
            fprintf(fp_bm, "%.2f %.6e %.6e\n", error_arr[baseline_idx], minsum_ber[baseline_idx], minsum_fer[baseline_idx]);
            fclose(fp_bm);
        }
        FILE *fp_bs = fopen(BASELINE_SPA_FILE, "w");
        if (fp_bs) {
            fprintf(fp_bs, "%.2f %.6e %.6e\n", error_arr[baseline_idx], spa_ber[baseline_idx], spa_fer[baseline_idx]);
            fclose(fp_bs);
        }
    }

    // ---- baseline 대비 배율 요약 (콘솔 + 파일) ----
    FILE *fp_sum = fopen(SUMMARY_FILE, "w");

    printf("\n========================================\n");
    printf(" Comparison vs Correct Estimation (error = 0)\n");
    printf("========================================\n");

    if (baseline_idx < 0) {
        printf("Warning: baseline(error=0) point not found in sweep range.\n");
    } else {
        double base_ms_ber = minsum_ber[baseline_idx], base_ms_fer = minsum_fer[baseline_idx];
        double base_spa_ber = spa_ber[baseline_idx], base_spa_fer = spa_fer[baseline_idx];

        printf("Error(dB)  MinSum_BER_x   MinSum_FER_x   SPA_BER_x      SPA_FER_x\n");
        if (fp_sum) {
            fprintf(fp_sum, "Error(dB) MinSum_BER MinSum_BER_ratio MinSum_FER MinSum_FER_ratio "
                            "SPA_BER SPA_BER_ratio SPA_FER SPA_FER_ratio\n");
        }

        for (int i = 0; i < count; i++) {
            double ms_ber_ratio = (base_ms_ber > 0.0) ? (minsum_ber[i] / base_ms_ber) : NAN;
            double ms_fer_ratio = (base_ms_fer > 0.0) ? (minsum_fer[i] / base_ms_fer) : NAN;
            double spa_ber_ratio = (base_spa_ber > 0.0) ? (spa_ber[i] / base_spa_ber) : NAN;
            double spa_fer_ratio = (base_spa_fer > 0.0) ? (spa_fer[i] / base_spa_fer) : NAN;

            printf("%+6.2f     %9.2fx     %9.2fx     %9.2fx     %9.2fx%s\n",
                   error_arr[i], ms_ber_ratio, ms_fer_ratio, spa_ber_ratio, spa_fer_ratio,
                   (i == baseline_idx) ? "   <- correct estimation (baseline)" : "");

            if (fp_sum) {
                fprintf(fp_sum, "%.2f %.6e %.4f %.6e %.4f %.6e %.4f %.6e %.4f\n",
                        error_arr[i],
                        minsum_ber[i], ms_ber_ratio, minsum_fer[i], ms_fer_ratio,
                        spa_ber[i], spa_ber_ratio, spa_fer[i], spa_fer_ratio);
            }
        }
    }

    if (fp_sum) {
        fclose(fp_sum);
        printf("\nSaved to %s\n", SUMMARY_FILE);
    }
}

// gnuplot으로 Min-Sum(왼쪽)과 Exact SPA(오른쪽)를 나란히 비교하는 1x2 multiplot.
// - baseline(error=0) 지점은 큰 별 마커(pt 3, ps 3)로 강조.
// - baseline BER 값에 수평 점선을 그어서, 다른 오차 지점이 거기서 얼마나
//   벗어나는지 한눈에 비교할 수 있게 한다.
// - 한글 라벨은 Windows gnuplot에서 인코딩이 깨질 수 있어 전부 영어로 표기한다.
static void run_gnuplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Warning: failed to launch gnuplot. Is it installed and in PATH?\n");
        return;
    }

    fprintf(gp, "set terminal windows size 1400, 650\n");
    fprintf(gp, "set multiplot layout 1,2 title "
                "'Polar Code (N=1024, K=512): Effect of sigma^2 Estimation Error (True SNR = %.1f dB)' font ',13'\n",
            TRUE_SNR_DB);
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'SNR Estimation Error (dB) = Assumed SNR - True SNR'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp, "set arrow from 0, graph 0 to 0, graph 1 nohead lc rgb 'gray' dt 2\n");
    fprintf(gp, "set label 'Correct Estimation (error=0)' at 0, graph 0.97 offset 1,0 tc rgb 'gray' font ',9'\n");

    // 두 데이터 파일(BER, FER)의 실제 최소/최대값을 각 패널별로 따로 읽어서
    // y축 범위를 데이터에 딱 맞게 좁힌다. 10의 거듭제곱 단위로 자동 스케일하면
    // 값이 좁은 범위에 몰려 있을 때 위아래로 빈 공간이 많이 남아 곡선의 굴곡이
    // 잘 안 보이기 때문. Min-Sum과 Exact SPA의 절대적인 오류율 수준이 서로
    // 다를 수 있으므로, 한쪽에 맞춘 범위를 공유하지 않고 패널마다 각자 최적의
    // 범위를 따로 계산해서 적용한다.
    fprintf(gp, "stats '" MINSUM_FILE "' using 2 nooutput prefix 'A'\n");
    fprintf(gp, "stats '" MINSUM_FILE "' using 3 nooutput prefix 'B'\n");
    fprintf(gp, "stats '" SPA_FILE "' using 2 nooutput prefix 'C'\n");
    fprintf(gp, "stats '" SPA_FILE "' using 3 nooutput prefix 'D'\n");
    fprintf(gp, "ymin_ms = (A_min < B_min) ? A_min : B_min\n");
    fprintf(gp, "ymax_ms = (A_max > B_max) ? A_max : B_max\n");
    fprintf(gp, "ymin_spa = (C_min < D_min) ? C_min : D_min\n");
    fprintf(gp, "ymax_spa = (C_max > D_max) ? C_max : D_max\n");

    // ---- 왼쪽 패널: Min-Sum ----
    fprintf(gp, "set yrange [ymin_ms*0.8:ymax_ms*1.3]\n");
    fprintf(gp, "set title 'Min-Sum (Linear -> Scale-Invariant, expected FLAT)'\n");
    fprintf(gp,
            "plot '" MINSUM_FILE "' using 1:2 with linespoints lw 2 pt 7 ps 0.8 lc rgb 'red' title 'BER', "
            "'" MINSUM_FILE "' using 1:3 with linespoints lw 2 pt 7 ps 0.8 lc rgb 'blue' title 'FER', "
            "'" BASELINE_MINSUM_FILE "' using 1:2 with points pt 3 ps 3 lc rgb 'red' title 'BER (baseline)', "
            "'" BASELINE_MINSUM_FILE "' using 1:3 with points pt 3 ps 3 lc rgb 'blue' title 'FER (baseline)'\n");

    // ---- 오른쪽 패널: Exact SPA ----
    fprintf(gp, "set yrange [ymin_spa*0.8:ymax_spa*1.3]\n");
    fprintf(gp, "set title 'Exact SPA (Nonlinear -> Sensitive to Error)'\n");
    fprintf(gp,
            "plot '" SPA_FILE "' using 1:2 with linespoints lw 2 pt 7 ps 0.8 lc rgb 'red' title 'BER', "
            "'" SPA_FILE "' using 1:3 with linespoints lw 2 pt 7 ps 0.8 lc rgb 'blue' title 'FER', "
            "'" BASELINE_SPA_FILE "' using 1:2 with points pt 3 ps 3 lc rgb 'red' title 'BER (baseline)', "
            "'" BASELINE_SPA_FILE "' using 1:3 with points pt 3 ps 3 lc rgb 'blue' title 'FER (baseline)'\n");

    fprintf(gp, "unset multiplot\n");
    PCLOSE(gp);
}

int main(void) {
    // 참고: rng_state는 run_monte_carlo() 안에서 매번 COMMON_SEED로 다시
    // 초기화되므로, 여기서 별도로 seed_xorshift64()를 호출할 필요가 없습니다
    // (공통 난수 기법을 쓰기 때문에 실행할 때마다 항상 같은 결과가 나옵니다).
    ensure_result_dir();

    printf("Running Polar Code (N=%d, K=%d) sigma^2 estimation error sweep "
           "(Min-Sum vs Exact SPA)...\n", N, K);
    run_error_sweep();

    printf("\nSimulation completed!\n");

    run_gnuplot();
    printf("Gnuplot window launched (if gnuplot is installed).\n");

    return 0;
}