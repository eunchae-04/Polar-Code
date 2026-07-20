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

// =================================================================
// 고성능 64비트 난수 생성기 (Xorshift64 Engine)
// =================================================================
uint64_t rng_state;

void seed_xorshift64() {
    rng_state = (uint64_t)time(NULL) ^ 0x5DEECE66DL;
    if (rng_state == 0) rng_state = 1;
}

uint64_t xorshift64() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

double rand_double() {
    return (double)xorshift64() / 18446744073709551615ULL;
}

// =================================================================
// Polar Code 구성 및 복호화 함수 (Exact SPA 적용)
// =================================================================
double phi(double x) {
    if (x <= 0.0) return 1.0;
    if (x <= 10.0) return exp(-0.4527 * pow(x, 0.86) + 0.0218);
    return sqrt(M_PI / x) * (1.0 - 10.0 / (7.0 * x)) * exp(-x / 4.0);
}

double phi_inv(double y) {
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

void construct_frozen_mask(double sigma2, int *info_mask) {
    double mu[N];
    mu[0] = 2.0 / sigma2;
    int n = 10;
    for (int stage = 1; stage <= n; stage++) {
        int block_size = 1 << (stage - 1);
        for (int i = 0; i < block_size; i++) {
            double T = mu[i];
            mu[i] = phi_inv(1.0 - pow(1.0 - phi(T), 2));
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
        int tmp = idx[i]; idx[i] = idx[max_j]; idx[max_j] = tmp;
    }
    for (int i = 0; i < N; i++) info_mask[i] = 0;
    for (int i = 0; i < K; i++) info_mask[idx[i]] = 1;
}

double gauss_box_muller() {
    double u1, u2;
    do { u1 = rand_double(); } while (u1 <= 1e-12);
    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

void polar_encode_recursive(const int *u, int *x, int n_len) {
    if (n_len == 1) {
        x[0] = u[0];
        return;
    }
    int n2 = n_len / 2;
    int v[n2];
    for (int i = 0; i < n2; i++) v[i] = u[i] ^ u[i + n2];
    polar_encode_recursive(v, x, n2);
    polar_encode_recursive(u + n2, x + n2, n2);
}

void polar_sc_decode_recursive(const double *llr, const int *info_mask, int *u_hat, int *u_coded, int n_len, int *mask_offset) {
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
        llr_left[i] = sign * min_val + log(1.0 + exp(-(a1 + a2))) - log(1.0 + exp(-fabs(a1 - a2)));
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
// 메인 시뮬레이션 컨트롤러
// =================================================================
int main(void) {
    seed_xorshift64();
    
    const double EbNo_start = 0.0;
    const double EbNo_end = 3.0;
    const double EbNo_step = 0.25;
    
    // -------------------------------------------------------------
    // [루프 1] Baseline 시나리오 (오차 없는 최적 환경)
    // -------------------------------------------------------------
    FILE *fp1 = fopen("simulation_baseline.txt", "w");
    printf("\n[1/3] Running Baseline Scenario...\n");
    printf("Eb/No (dB)\tBER\t\tFER\n");
    for (double snr_db = EbNo_start; snr_db <= EbNo_end + 1e-9; snr_db += EbNo_step) {
        double EbNo_linear = pow(10.0, snr_db / 10.0);
        double sigma2_true = 1.0 / (2.0 * R * EbNo_linear);
        
        int info_mask[N];
        construct_frozen_mask(sigma2_true, info_mask); // 채널에 최적화된 마스크 생성
        
        long total_errors = 0, total_bits = 0, total_frames = 0, total_frame_errors = 0;
        while (total_errors < 500) { // 연산 속도를 확보하기 위해 타겟 에러를 500개로 조정
            int info_bits[K];
            for (int i = 0; i < K; i++) info_bits[i] = (int)(xorshift64() & 1);
            int u[N], info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) u[i] = info_bits[info_ptr++];
                else u[i] = 0;
            }
            int x[N];
            polar_encode_recursive(u, x, N);
            double tx[N], rx[N];
            for (int i = 0; i < N; i++) tx[i] = 1.0 - 2.0 * x[i];
            for (int i = 0; i < N; i++) rx[i] = tx[i] + gauss_box_muller() * sqrt(sigma2_true);
            
            // LLR 오차 없음 (sigma2_true 사용)
            double llr[N];
            for (int i = 0; i < N; i++) llr[i] = (2.0 / sigma2_true) * rx[i];
            
            int u_hat[N], u_coded[N], mask_offset = 0;
            polar_sc_decode_recursive(llr, info_mask, u_hat, u_coded, N, &mask_offset);
            int bit_errors = 0; info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) {
                    if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
                }
            }
            if (bit_errors > 0) total_frame_errors++;
            total_errors += bit_errors; total_bits += K; total_frames++;
        }
        double final_ber = (double)total_errors / (double)total_bits;
        double final_fer = (double)total_frame_errors / (double)total_frames;
        printf("%.2f dB\t\t%.6e\t%.6e\n", snr_db, final_ber, final_fer);
        fprintf(fp1, "%.2f %.6e %.6e\n", snr_db, final_ber, final_fer);
    }
    fclose(fp1);

    // -------------------------------------------------------------
    // [루프 2] 시나리오 A (LLR 스케일링 오차, alpha = 0.5)
    // -------------------------------------------------------------
    FILE *fp2 = fopen("simulation_scenario_a.txt", "w");
    double alpha = 0.5; // 잡음을 반으로 과소평가하여 LLR 스케일이 비정상적으로 커진 상태
    printf("\n[2/3] Running Scenario A (LLR Mismatch, alpha = 0.5)...\n");
    printf("Eb/No (dB)\tBER\t\tFER\n");
    for (double snr_db = EbNo_start; snr_db <= EbNo_end + 1e-9; snr_db += EbNo_step) {
        double EbNo_linear = pow(10.0, snr_db / 10.0);
        double sigma2_true = 1.0 / (2.0 * R * EbNo_linear);
        double sigma2_est = alpha * sigma2_true; // 수신기가 오추정한 잡음 분산
        
        int info_mask[N];
        construct_frozen_mask(sigma2_true, info_mask);
        
        long total_errors = 0, total_bits = 0, total_frames = 0, total_frame_errors = 0;
        while (total_errors < 500) {
            int info_bits[K];
            for (int i = 0; i < K; i++) info_bits[i] = (int)(xorshift64() & 1);
            int u[N], info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) u[i] = info_bits[info_ptr++];
                else u[i] = 0;
            }
            int x[N];
            polar_encode_recursive(u, x, N);
            double tx[N], rx[N];
            for (int i = 0; i < N; i++) tx[i] = 1.0 - 2.0 * x[i];
            for (int i = 0; i < N; i++) rx[i] = tx[i] + gauss_box_muller() * sqrt(sigma2_true);
            
            // LLR 계산 시 왜곡된 잡음 분산(sigma2_est) 사용으로 Jacobian 보정항 교란 유도
            double llr[N];
            for (int i = 0; i < N; i++) llr[i] = (2.0 / sigma2_est) * rx[i];
            
            int u_hat[N], u_coded[N], mask_offset = 0;
            polar_sc_decode_recursive(llr, info_mask, u_hat, u_coded, N, &mask_offset);
            int bit_errors = 0; info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) {
                    if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
                }
            }
            if (bit_errors > 0) total_frame_errors++;
            total_errors += bit_errors; total_bits += K; total_frames++;
        }
        double final_ber = (double)total_errors / (double)total_bits;
        double final_fer = (double)total_frame_errors / (double)total_frames;
        printf("%.2f dB\t\t%.6e\t%.6e\n", snr_db, final_ber, final_fer);
        fprintf(fp2, "%.2f %.6e %.6e\n", snr_db, final_ber, final_fer);
    }
    fclose(fp2);

    // -------------------------------------------------------------
    // [루프 3] 시나리오 B (Design SNR 미스매치, Design SNR = 0.0 dB 고정)
    // -------------------------------------------------------------
    FILE *fp3 = fopen("simulation_scenario_b.txt", "w");
    double design_snr_db = 0.0; // 부호 설계용 SNR을 0.0 dB로 완전히 박아둠
    double design_EbNo_linear = pow(10.0, design_snr_db / 10.0);
    double sigma2_design = 1.0 / (2.0 * R * design_EbNo_linear);
    int fixed_info_mask[N];
    construct_frozen_mask(sigma2_design, fixed_info_mask); // ★ 루프 밖에서 0dB 전용 마스크 딱 하나만 생성
    
    printf("\n[3/3] Running Scenario B (Design SNR Mismatch = 0.0 dB Fixed)...\n");
    printf("Eb/No (dB)\tBER\t\tFER\n");
    for (double snr_db = EbNo_start; snr_db <= EbNo_end + 1e-9; snr_db += EbNo_step) {
        double EbNo_linear = pow(10.0, snr_db / 10.0);
        double sigma2_true = 1.0 / (2.0 * R * EbNo_linear);
        
        long total_errors = 0, total_bits = 0, total_frames = 0, total_frame_errors = 0;
        while (total_errors < 500) {
            int info_bits[K];
            for (int i = 0; i < K; i++) info_bits[i] = (int)(xorshift64() & 1);
            int u[N], info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (fixed_info_mask[i]) u[i] = info_bits[info_ptr++]; // ★ 고정된 0dB용 마스크 적용
                else u[i] = 0;
            }
            int x[N];
            polar_encode_recursive(u, x, N);
            double tx[N], rx[N];
            for (int i = 0; i < N; i++) tx[i] = 1.0 - 2.0 * x[i];
            for (int i = 0; i < N; i++) rx[i] = tx[i] + gauss_box_muller() * sqrt(sigma2_true);
            
            double llr[N];
            for (int i = 0; i < N; i++) llr[i] = (2.0 / sigma2_true) * rx[i];
            
            int u_hat[N], u_coded[N], mask_offset = 0;
            polar_sc_decode_recursive(llr, fixed_info_mask, u_hat, u_coded, N, &mask_offset);
            int bit_errors = 0; info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (fixed_info_mask[i]) {
                    if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
                }
            }
            if (bit_errors > 0) total_frame_errors++;
            total_errors += bit_errors; total_bits += K; total_frames++;
        }
        double final_ber = (double)total_errors / (double)total_bits;
        double final_fer = (double)total_frame_errors / (double)total_frames;
        printf("%.2f dB\t\t%.6e\t%.6e\n", snr_db, final_ber, final_fer);
        fprintf(fp3, "%.2f %.6e %.6e\n", snr_db, final_ber, final_fer);
    }
    fclose(fp3);

    // =============================================================
    // [Gnuplot 연동] 화면 분할(Multiplot) 동적 그래프 시각화 명령
    // =============================================================
    FILE *gp = popen("gnuplot -persist", "w");
    if (gp != NULL) {
        // 창 크기를 가로로 길게 배치 (3개 플롯을 늘어놓기 위함)
        fprintf(gp, "set terminal windows size 1300, 450\n"); 
        
        // 1행 3열 구조의 멀티플롯 레이아웃 선언
        fprintf(gp, "set multiplot layout 1,3 title 'Polar Code Scenario Analyses (N=1024, K=512)' font ',13'\n");
        
        // 공통 속성 정의
        fprintf(gp, "set logscale y\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set yrange [1e-4:1]\n");
        fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
        fprintf(gp, "set ylabel 'Error Probability'\n");
        
        // 1번째 서브플롯: Baseline
        fprintf(gp, "set title '1. Baseline (Perfect)'\n");
        fprintf(gp, "plot 'simulation_baseline.txt' using 1:2 with linespoints lw 2 lc rgb 'purple' title 'BER', 'simulation_baseline.txt' using 1:3 with linespoints lw 2 lc rgb 'cyan' title 'FER'\n");
        
        // 2번째 서브플롯: Scenario A (Y축 라벨 숨겨서 깔끔하게 정리)
        fprintf(gp, "set title '2. Scenario A (LLR Error, a=0.5)'\n");
        fprintf(gp, "plot 'simulation_scenario_a.txt' using 1:2 with linespoints lw 2 lc rgb 'red' title 'BER', 'simulation_scenario_a.txt' using 1:3 with linespoints lw 2 lc rgb 'orange' title 'FER'\n");
        
        // 3번째 서브플롯: Scenario B
        fprintf(gp, "set title '3. Scenario B (Design SNR=0dB)'\n");
        fprintf(gp, "plot 'simulation_scenario_b.txt' using 1:2 with linespoints lw 2 lc rgb 'blue' title 'BER', 'simulation_scenario_b.txt' using 1:3 with linespoints lw 2 lc rgb 'dark-green' title 'FER'\n");
        
        // 멀티플롯 종료 선언
        fprintf(gp, "unset multiplot\n");
        pclose(gp);
    }

    printf("\nAll simulations completed! Gnuplot multiplot window popped up.\n");
    return 0;
}