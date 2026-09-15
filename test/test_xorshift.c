/*
 * Xorshift64 균일 분포 검증 스크립트
 * 컴파일: gcc -O2 -o test_xorshift test_xorshift.c -lm
 * 실행: ./test_xorshift
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#define NUM_SAMPLES 100000000 // 1억 회 추출
#define NUM_BINS 20           // [0.0, 1.0) 구간을 20등분 (간격: 0.05)
#define DATA_FILE "xorshift_dist.txt"

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

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

static void plot_histogram(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Notice: Gnuplot이 설치되어 있지 않거나 실행할 수 없어 그래프 생성을 건너뜁니다.\n");
        return;
    }

    fprintf(gp, "set title 'Xorshift64 Uniform Distribution Test (100 Million Samples)'\n");
    fprintf(gp, "set xlabel 'Value Range [0.0, 1.0)'\n");
    fprintf(gp, "set ylabel 'Sample Count'\n");
    fprintf(gp, "set grid y\n");
    fprintf(gp, "set style fill solid 0.6 border -1\n");
    fprintf(gp, "set boxwidth 0.04\n");
    fprintf(gp, "plot '%s' using 1:2 with boxes lc rgb 'steelblue' title 'Bin Count'\n", DATA_FILE);

    PCLOSE(gp);
}

int main(void) {
    seed_xorshift64();

    uint64_t bin_counts[NUM_BINS] = {0};
    double sum = 0.0;
    double sum_sq = 0.0;

    printf("1억 개의 샘플을 생성하여 Xorshift64 분포를 검증 중입니다...\n");

    for (uint64_t i = 0; i < NUM_SAMPLES; i++) {
        double val = rand_double();

        // 수치 안정성 경계 처리
        if (val >= 1.0) val = 0.999999999;
        if (val < 0.0) val = 0.0;

        sum += val;
        sum_sq += val * val;

        int bin_idx = (int)(val * NUM_BINS);
        if (bin_idx >= NUM_BINS) bin_idx = NUM_BINS - 1;
        bin_counts[bin_idx]++;
    }

    // 1. 통계치 측정
    double sample_mean = sum / (double)NUM_SAMPLES;
    double sample_variance = (sum_sq / (double)NUM_SAMPLES) - (sample_mean * sample_mean);

    // 이론값: [0, 1] 균일분포 평균 = 0.5, 분산 = 1/12 ≈ 0.083333333
    double theoretical_mean = 0.5;
    double theoretical_variance = 1.0 / 12.0;

    printf("\n=======================================================\n");
    printf("                  통계적 검증 결과                     \n");
    printf("=======================================================\n");
    printf("총 샘플 수       : %u 회\n", NUM_SAMPLES);
    printf("측정 평균 (Mean) : %.8f (이론값: %.8f, 오차: %.8f)\n", 
           sample_mean, theoretical_mean, fabs(sample_mean - theoretical_mean));
    printf("측정 분산 (Var)  : %.8f (이론값: %.8f, 오차: %.8f)\n", 
           sample_variance, theoretical_variance, fabs(sample_variance - theoretical_variance));
    printf("=======================================================\n\n");

    // 2. 파일 저장
    FILE *fp = fopen(DATA_FILE, "w");
    if (fp == NULL) {
        printf("데이터 파일 생성 실패\n");
        return 1;
    }

    double expected_per_bin = (double)NUM_SAMPLES / (double)NUM_BINS;
    printf("구간별 빈도수 (이상적인 구간당 샘플 수: %.0f 개)\n", expected_per_bin);
    printf("-------------------------------------------------------\n");
    printf("  구간 (Bin Range)      |   샘플 수 (Count)   | 비율 (%)\n");
    printf("-------------------------------------------------------\n");

    for (int i = 0; i < NUM_BINS; i++) {
        double bin_center = (i + 0.5) / (double)NUM_BINS;
        double ratio = ((double)bin_counts[i] / (double)NUM_SAMPLES) * 100.0;
        
        printf("  [%.2f ~ %.2f)       |   %10lu    |  %.3f%%\n", 
               (double)i / NUM_BINS, (double)(i + 1) / NUM_BINS, bin_counts[i], ratio);
        
        fprintf(fp, "%.4f %lu\n", bin_center, bin_counts[i]);
    }
    printf("-------------------------------------------------------\n");
    fclose(fp);

    // 3. Gnuplot 시각화
    plot_histogram();

    return 0;
}