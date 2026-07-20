#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

// ★ [쥬피터 환경 대응 패치] M_PI 가 선언되어 있지 않은 경우 직접 매핑하여 에러 원천 차단
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 1024
#define K 512
#define R ((double)K / (double)N)

// [GA 구축] 정석 phi 함수 구현
double phi(double x) {
    if (x <= 0.0) return 1.0;
    if (x <= 10.0) return exp(-0.4527 * pow(x, 0.86) + 0.0218);
    return sqrt(M_PI / x) * (1.0 - 10.0 / (7.0 * x)) * exp(-x / 4.0);
}

// [GA 구축] 상한선/하한선 제한이 해제된 정밀 이진 탐색 phi_inv
double phi_inv(double y) {
    if (y >= 1.0) return 0.0;
    if (y < 1e-200) y = 1e-200; // 언더플로우 최소 한계 보정
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

// [GA 구축] 각 SNR 단계별 정밀 Frozen Mask 생성 함수
void construct_frozen_mask(double sigma2, int *info_mask) {
    double mu[N];
    mu[0] = 2.0 / sigma2;
    
    int n = 10; // log2(1024)
    for (int stage = 1; stage <= n; stage++) {
        int block_size = 1 << (stage - 1);
        for (int i = 0; i < block_size; i++) {
            double T = mu[i];
            mu[i] = phi_inv(1.0 - pow(1.0 - phi(T), 2)); // Bad Channel (W-)
            mu[i + block_size] = 2.0 * T;                // Good Channel (W+)
        }
    }
    
    // Bottom-Up으로 계산된 신뢰도를 Natural Order 물리 트리로 매핑 (Bit-Reversal)
    double natural_mu[N];
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < n; j++) {
            if (i & (1 << j)) rev |= (1 << (n - 1 - j));
        }
        natural_mu[rev] = mu[i];
    }
    
    // 신뢰도 내림차순 인덱스 정렬
    int idx[N];
    for (int i = 0; i < N; i++) idx[i] = i;
    for (int i = 0; i < N - 1; i++) {
        int max_j = i;
        for (int j = i + 1; j < N; j++) {
            if (natural_mu[idx[j]] > natural_mu[idx[max_j]]) max_j = j;
        }
        int tmp = idx[i]; idx[i] = idx[max_j]; idx[max_j] = tmp;
    }
    
    // 상위 K개는 Information(1), 나머지는 Frozen(0)
    for (int i = 0; i < N; i++) info_mask[i] = 0;
    for (int i = 0; i < K; i++) info_mask[idx[i]] = 1;
}

// [채널] Box-Muller 변환 기반 AWGN 생성기
double gauss_box_muller() {
    double u1, u2;
    do { u1 = (double)rand() / (double)RAND_MAX; } while (u1 <= 1e-12);
    u2 = (double)rand() / (double)RAND_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// [송신단] Top-Down 재귀형 Natural Order 인코더
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

// [수신단] 무결성이 검증된 진짜 1024-SC 재귀 복호화기 (Plain Concat 조합)
void polar_sc_decode_recursive(const double *llr, const int *info_mask, int *u_hat, int *u_coded, int n_len, int *mask_offset) {
    if (n_len == 1) {
        if (info_mask[*mask_offset]) {
            u_hat[0] = (llr[0] < 0.0) ? 1 : 0; // Information 비트 경판정
        } else {
            u_hat[0] = 0;                      // Frozen 비트 강제 영 고정
        }
        u_coded[0] = u_hat[0];
        (*mask_offset)++;
        return;
    }
    
    int n2 = n_len / 2;
    double llr_left[n2];
    double llr_right[n2];
    int u_hat_left[n2], u_hat_right[n2];
    int u_coded_left[n2], u_coded_right[n2];
    
    // 1. Left Path: f-function 연산 (Min-Sum)
    for (int i = 0; i < n2; i++) {
        double l1 = llr[i];
        double l2 = llr[i + n2];
        double sign = ((l1 < 0.0) ^ (l2 < 0.0)) ? -1.0 : 1.0;
        double a1 = fabs(l1);
        double a2 = fabs(l2);
        llr_left[i] = sign * (a1 < a2 ? a1 : a2);
    }
    
    polar_sc_decode_recursive(llr_left, info_mask, u_hat_left, u_coded_left, n2, mask_offset);
    
    // 2. Right Path: g-function 연산 (결정 피드백 반영)
    for (int i = 0; i < n2; i++) {
        llr_right[i] = llr[i + n2] + (1 - 2 * u_coded_left[i]) * llr[i];
    }
    
    polar_sc_decode_recursive(llr_right, info_mask, u_hat_right, u_coded_right, n2, mask_offset);
    
    // 3. 상위 계층 코드로의 병합 및 반환 (Plain Concat 매핑)
    for (int i = 0; i < n2; i++) {
        // 정보 비트(u_hat)는 순수 단순 연결로 결합합니다.
        u_hat[i] = u_hat_left[i];                  
        u_hat[i + n2] = u_hat_right[i];
        
        // 채널 Coded 비트에 인코더 대칭 상위 기하학 구조(XOR)를 반영합니다.
        u_coded[i] = u_coded_left[i] ^ u_coded_right[i]; 
        u_coded[i + n2] = u_coded_right[i];
    }   
}

// [메인 컨트롤러] 몬테카를로 루프 엔진
int main(void) {
    srand((unsigned)time(NULL));
    
    const double EbNo_start = 0.0;
    const double EbNo_end = 3.0;
    const double EbNo_step = 0.25;
    
    // 데이터를 저장할 텍스트 파일 열기
    FILE *fp = fopen("simulation_results.txt", "w");
    if (fp == NULL) {
        printf("오류: 데이터 파일을 생성할 수 없습니다.\n");
        return -1;
    }
    
    printf("=============================================\n");
    printf("Eb/No (dB)\tBER\t\tFER\n");
    printf("=============================================\n");
    
    for (double snr_db = EbNo_start; snr_db <= EbNo_end + 1e-9; snr_db += EbNo_step) {
        double EbNo_linear = pow(10.0, snr_db / 10.0);
        double sigma2 = 1.0 / (2.0 * R * EbNo_linear);
        
        // SNR 분산에 대응하는 정확한 가우시안 근사 마스크 형성
        int info_mask[N];
        construct_frozen_mask(sigma2, info_mask);
        
        long total_errors = 0;
        long total_bits = 0;
        long total_frames = 0;
        long total_frame_errors = 0;
        
        while (total_errors < 1000) {
            // 1) 무작위 정보 비트 K개 생성
            int info_bits[K];
            for (int i = 0; i < K; i++) info_bits[i] = rand() & 1;
            
            // 2) 마스크 채널 매핑 순서에 입각하여 u 벡터 조립
            int u[N];
            int info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) u[i] = info_bits[info_ptr++];
                else u[i] = 0;
            }
            
            // 3) 인코더 가동
            int x[N];
            polar_encode_recursive(u, x, N);
            
            // 4) BPSK 변조 및 AWGN 채널 통과
            double tx[N], rx[N];
            for (int i = 0; i < N; i++) tx[i] = 1.0 - 2.0 * x[i];
            for (int i = 0; i < N; i++) {
                rx[i] = tx[i] + gauss_box_muller() * sqrt(sigma2);
            }
            
            // 5) 수신단 초기 LLR 도출
            double llr[N];
            for (int i = 0; i < N; i++) llr[i] = (2.0 / sigma2) * rx[i];
            
            // 6) 무결성 재귀 디코더 연산
            int u_hat[N], u_coded[N];
            int mask_offset = 0;
            polar_sc_decode_recursive(llr, info_mask, u_hat, u_coded, N, &mask_offset);
            
            // 7) 몬테카를로 통계치 누적 집계
            int bit_errors = 0;
            info_ptr = 0;
            for (int i = 0; i < N; i++) {
                if (info_mask[i]) {
                    if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
                }
            }
            
            if (bit_errors > 0) total_frame_errors++;
            total_errors += bit_errors;
            total_bits += K;
            total_frames++;
        }
        
        double final_ber = (double)total_errors / (double)total_bits;
        double final_fer = (double)total_frame_errors / (double)total_frames;
        
        // 콘솔창에 결과 출력
        printf("%.2f dB\t\t%.6e\t%.6e\n", snr_db, final_ber, final_fer);
        fflush(stdout);
        
        // 생성한 텍스트 파일에 실시간 데이터 입력 (공백 기준 분리)
        fprintf(fp, "%.2f %.6e %.6e\n", snr_db, final_ber, final_fer);
    }
    printf("=============================================\n");

    // 시뮬레이션 종료 후 데이터 파일 닫기
    fclose(fp);

    // Gnuplot을 실행하여 방금 저장한 텍스트 파일 불러와서 그리기
    FILE *gp = popen("gnuplot -persist", "w");
    if (gp != NULL) {
        fprintf(gp, "set title 'Polar Code (N=1024, K=512) SC Performance'\n");
        fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
        fprintf(gp, "set ylabel 'Error Probability'\n");
        fprintf(gp, "set logscale y\n"); // Y축 로그 스케일 설정
        fprintf(gp, "set grid\n");
        
        // 하드코딩 대신 파일(simulation_results.txt)을 직접 지정하여 플롯 명령 수행
        // 1열(Eb/No) 기준 2열(BER), 1열 기준 3열(FER)을 선과 점(linespoints)으로 매핑
        fprintf(gp, "plot 'simulation_results.txt' using 1:2 with linespoints lw 2 title 'BER', 'simulation_results.txt' using 1:3 with linespoints lw 2 title 'FER'\n");
        
        pclose(gp);
    }

    return 0;
}