#ifndef POLAR_CONFIG_H
#define POLAR_CONFIG_H

#include <stdbool.h>

#define N 1024
#define K 512
#define R ((double)K / (double)N)

#define EBNO_START 0.0
#define EBNO_END   3.0
#define EBNO_STEP  0.25
#define TARGET_ERRORS 10000

typedef struct {
    const char *output_file;
    const char *title;
    double llr_scale;        // sigma2_est = sigma2_true * llr_scale
    bool use_fixed_mask;      // fixed_mask_snr_db로 생성한 마스크를 고정 사용
    double fixed_mask_snr_db; // fixed mask 생성용 SNR
} SimulationConfig;

#endif