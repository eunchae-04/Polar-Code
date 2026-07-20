#include "polar_math.h"
#include "polar_config.h"

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

double compute_sigma2_from_ebno_db(double ebno_db) {
    double ebno_linear = pow(10.0, ebno_db / 10.0);
    return 1.0 / (2.0 * R * ebno_linear);
}

void construct_frozen_mask(double sigma2, int *info_mask) {
    double mu[N];
    mu[0] = 2.0 / sigma2;

    const int n = 10; // N = 2^10 = 1024
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
