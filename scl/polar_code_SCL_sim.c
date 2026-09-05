/*
    gcc -O2 -Wall -o polar_code_SCL_sim polar_code_SCL_sim.c -lm
    ./polar_code_SCL_sim
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

// SCL 리스트 크기 (L). 크게 할수록 성능은 좋아지지만 계산량은 대략 L배로 늘어납니다.
// 보통 2, 4, 8, 16, 32 중에서 선택합니다.
#define LIST_SIZE 8

#define EBNO_START 0.0
#define EBNO_END   3.0
#define EBNO_STEP  0.25
#define TARGET_ERRORS 10000

#define RESULT_DIR "result"
#define RESULT_FILE RESULT_DIR "/simulation_results_scl.txt"

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
// Encode
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

// =================================================================
// SCL (Successive Cancellation List) 디코더
// =================================================================
// 각 경로(path)의 누적 Path Metric. 값이 작을수록 신뢰도가 높은 경로.
// PM 갱신 규칙(간이 근사식, Balatsoukas-Stimming 등에서 널리 쓰이는 형태):
//   - 결정한 비트가 LLR 부호가 가리키는 하드 디시전과 일치하면 페널티 없음
//   - 불일치하면 |LLR| 만큼 페널티 추가
static double path_metric[LIST_SIZE];

// n_len 단위 재귀 함수. 진입 시 *list_size 개의 경로가 활성 상태이며,
// leaf(정보 비트)에서 최대 2배로 늘어났다가 LIST_SIZE개로 가지치기(pruning)된다.
// origin[q]는 "새 경로 q가 이 함수에 들어올 때의 어느 경로(0..list_size_in-1)에서
// 갈라져 나왔는지"를 상위(부모) 재귀 호출이 참조할 수 있도록 반환하는 값이다.
static void scl_decode_recursive(int n_len,
                                  double llr[LIST_SIZE][n_len],
                                  const int *info_mask,
                                  int u_hat[LIST_SIZE][n_len],
                                  int u_coded[LIST_SIZE][n_len],
                                  int *mask_offset,
                                  int *list_size,
                                  int origin[LIST_SIZE]) {
    if (n_len == 1) {
        int in_size = *list_size;

        if (info_mask[*mask_offset]) {
            // 정보 비트: 각 경로를 0/1 두 가지로 분기시켜 최대 2*in_size개의
            // 후보를 만든 뒤, Path Metric이 작은 순서로 LIST_SIZE개만 남긴다.
            int cand_parent[2 * LIST_SIZE];
            int cand_bit[2 * LIST_SIZE];
            double cand_metric[2 * LIST_SIZE];
            int cand_count = 0;

            for (int p = 0; p < in_size; p++) {
                double l = llr[p][0];
                double pm = path_metric[p];

                cand_parent[cand_count] = p;
                cand_bit[cand_count] = 0;
                cand_metric[cand_count] = pm + (l < 0.0 ? fabs(l) : 0.0);
                cand_count++;

                cand_parent[cand_count] = p;
                cand_bit[cand_count] = 1;
                cand_metric[cand_count] = pm + (l >= 0.0 ? fabs(l) : 0.0);
                cand_count++;
            }

            // Path Metric 오름차순 정렬 (후보 수가 최대 2*LIST_SIZE개라 선택 정렬로 충분)
            int order[2 * LIST_SIZE];
            for (int i = 0; i < cand_count; i++) order[i] = i;
            for (int i = 0; i < cand_count - 1; i++) {
                int min_j = i;
                for (int j = i + 1; j < cand_count; j++) {
                    if (cand_metric[order[j]] < cand_metric[order[min_j]]) min_j = j;
                }
                int tmp = order[i];
                order[i] = order[min_j];
                order[min_j] = tmp;
            }

            int new_size = (cand_count < LIST_SIZE) ? cand_count : LIST_SIZE;

            // 기존 배열을 바로 덮어쓰면 같은 부모가 두 자식 모두에 선택될 때
            // 값이 꼬일 수 있으므로 임시 버퍼에 모은 뒤 한 번에 반영한다.
            int new_origin[LIST_SIZE];
            int new_bit[LIST_SIZE];
            double new_pm[LIST_SIZE];

            for (int q = 0; q < new_size; q++) {
                int c = order[q];
                new_origin[q] = cand_parent[c];
                new_bit[q] = cand_bit[c];
                new_pm[q] = cand_metric[c];
            }

            for (int q = 0; q < new_size; q++) {
                u_hat[q][0] = new_bit[q];
                u_coded[q][0] = new_bit[q];
                origin[q] = new_origin[q];
                path_metric[q] = new_pm[q];
            }

            *list_size = new_size;
        } else {
            // 프로즌 비트: 값은 무조건 0으로 고정, 경로 수는 그대로.
            // LLR 부호가 1을 가리키는데 억지로 0을 선택한 경우 페널티 반영.
            for (int p = 0; p < in_size; p++) {
                double l = llr[p][0];
                u_hat[p][0] = 0;
                u_coded[p][0] = 0;
                origin[p] = p;
                if (l < 0.0) path_metric[p] += fabs(l);
            }
        }

        (*mask_offset)++;
        return;
    }

    int n2 = n_len / 2;
    int list_size_in = *list_size;

    // f-function (Min-Sum) : 왼쪽 절반 LLR 계산, 경로별로 동일하게 적용
    double llr_left[LIST_SIZE][n2];
    for (int p = 0; p < list_size_in; p++) {
        for (int i = 0; i < n2; i++) {
            double l1 = llr[p][i];
            double l2 = llr[p][i + n2];
            double sign = ((l1 < 0.0) ^ (l2 < 0.0)) ? -1.0 : 1.0;
            double a1 = fabs(l1);
            double a2 = fabs(l2);
            double min_val = (a1 < a2) ? a1 : a2;
            llr_left[p][i] = sign * min_val;
        }
    }

    int u_hat_left[LIST_SIZE][n2];
    int u_coded_left[LIST_SIZE][n2];
    int origin_left[LIST_SIZE];

    scl_decode_recursive(n2, llr_left, info_mask, u_hat_left, u_coded_left,
                          mask_offset, list_size, origin_left);
    int list_size_mid = *list_size;

    // g-function (결정 피드백) : 왼쪽에서 갈라진 각 경로에 대해 오른쪽 LLR 계산.
    // llr[origin_left[q]]는 "이 경로가 이 함수에 들어올 때 갖고 있던" 원본 LLR.
    double llr_right[LIST_SIZE][n2];
    for (int q = 0; q < list_size_mid; q++) {
        int src = origin_left[q];
        for (int i = 0; i < n2; i++) {
            llr_right[q][i] = llr[src][n2 + i] + (1 - 2 * u_coded_left[q][i]) * llr[src][i];
        }
    }

    int u_hat_right[LIST_SIZE][n2];
    int u_coded_right[LIST_SIZE][n2];
    int origin_right[LIST_SIZE];

    scl_decode_recursive(n2, llr_right, info_mask, u_hat_right, u_coded_right,
                          mask_offset, list_size, origin_right);
    int list_size_final = *list_size;

    // 왼쪽/오른쪽 결과를 합쳐서 이 레벨의 u_hat, u_coded, origin을 완성.
    for (int q = 0; q < list_size_final; q++) {
        int r = origin_right[q]; // list_size_mid 공간에서의 인덱스
        origin[q] = origin_left[r]; // list_size_in 공간으로 합성(compose)
        for (int i = 0; i < n2; i++) {
            u_hat[q][i] = u_hat_left[r][i];
            u_hat[q][n2 + i] = u_hat_right[q][i];
            u_coded[q][i] = u_coded_left[r][i] ^ u_coded_right[q][i];
            u_coded[q][n2 + i] = u_coded_right[q][i];
        }
    }
}

// 채널 LLR(N개)을 받아 SCL 디코딩을 수행하고, Path Metric이 가장 작은 경로를
// 최종 결과(best_u_hat, N개, natural order)로 반환한다. (CRC 미사용, PM 최소 경로 선택)
static void polar_scl_decode(const double *channel_llr, const int *info_mask, int *best_u_hat) {
    static double llr0[LIST_SIZE][N];
    static int u_hat[LIST_SIZE][N];
    static int u_coded[LIST_SIZE][N];
    int origin[LIST_SIZE];

    for (int i = 0; i < N; i++) llr0[0][i] = channel_llr[i];
    path_metric[0] = 0.0;

    int list_size = 1;
    int mask_offset = 0;

    scl_decode_recursive(N, llr0, info_mask, u_hat, u_coded, &mask_offset, &list_size, origin);

    int best_q = 0;
    for (int q = 1; q < list_size; q++) {
        if (path_metric[q] < path_metric[best_q]) best_q = q;
    }

    for (int i = 0; i < N; i++) best_u_hat[i] = u_hat[best_q][i];
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

// BPSK 변조 + AWGN 채널
static void awgn_channel(const int *x, double sigma2_true, double *rx) {
    double sigma = sqrt(sigma2_true);
    for (int i = 0; i < N; i++) {
        double tx = 1.0 - 2.0 * x[i]; // 0 -> +1, 1 -> -1
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
// 시뮬레이션 (Eb/No 스윕 + 몬테카를로, BER/FER 동시 산출)
// =================================================================
static void run_simulation(FILE *fp) {
    printf("\nEb/No (dB)\tBER\t\tFER\n");

    for (double snr_db = EBNO_START; snr_db <= EBNO_END + 1e-9; snr_db += EBNO_STEP) {
        double sigma2_true = compute_sigma2_from_ebno_db(snr_db);
        double sigma2_est = sigma2_true; // 완벽한 SNR 추정 가정

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

            generate_info_bits(info_bits);
            build_u_from_mask(info_mask, info_bits, u);
            polar_encode_recursive(u, x, N);
            awgn_channel(x, sigma2_true, rx);
            compute_llr(rx, sigma2_est, llr);
            polar_scl_decode(llr, info_mask, u_hat);

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

// gnuplot을 통해 BER/FER Waterfall 곡선을 자동으로 표시
static void run_gnuplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Warning: failed to launch gnuplot. Is it installed and in PATH?\n");
        return;
    }

    fprintf(gp, "set terminal windows size 800, 600\n");
    fprintf(gp, "set title 'Polar Code (N=1024, K=512) SCL Decoder Performance (L=%d)'\n", LIST_SIZE);
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
    seed_xorshift64();
    ensure_result_dir();

    FILE *fp = fopen(RESULT_FILE, "w");
    if (fp == NULL) {
        printf("Failed to open output file: %s\n", RESULT_FILE);
        return 1;
    }

    printf("Running Polar Code (N=%d, K=%d) SCL simulation (List size L=%d)...\n", N, K, LIST_SIZE);
    run_simulation(fp);
    fclose(fp);

    printf("\nSimulation completed! Results saved to %s\n", RESULT_FILE);

    run_gnuplot();
    printf("Gnuplot window launched (if gnuplot is installed).\n");

    return 0;
}