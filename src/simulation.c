#include "simulation.h"

#include "polar_math.h"
#include "polar_codec.h"
#include "rng.h"

#include <stdio.h>
#include <math.h>

void generate_info_bits(int *info_bits) {
    for (int i = 0; i < K; i++) {
        info_bits[i] = (int)(xorshift64() & 1ULL);
    }
}

void build_u_from_mask(const int *mask, const int *info_bits, int *u) {
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

void compute_llr(const double *rx, double sigma2_est, double *llr) {
    for (int i = 0; i < N; i++) {
        llr[i] = (2.0 / sigma2_est) * rx[i];
    }
}

long count_bit_errors(const int *info_mask, const int *u_hat, const int *info_bits) {
    int info_ptr = 0;
    long bit_errors = 0;

    for (int i = 0; i < N; i++) {
        if (info_mask[i]) {
            if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
        }
    }
    return bit_errors;
}

void run_simulation(const SimulationConfig *cfg, const int *fixed_mask) {
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