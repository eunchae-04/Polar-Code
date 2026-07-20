#include "polar_codec.h"

#include <math.h>

static double log1pexp_neg(double x) {
    if (x > 50.0) return exp(-x);
    return log1p(exp(-x));
}

void polar_encode_recursive(const int *u, int *x, int n_len) {
    if (n_len == 1) {
        x[0] = u[0];
        return;
    }

    int n2 = n_len / 2;
    int v[n2];

    for (int i = 0; i < n2; i++) {
        v[i] = u[i] ^ u[i + n2];
    }

    polar_encode_recursive(v, x, n2);
    polar_encode_recursive(u + n2, x + n2, n2);
}

void polar_sc_decode_recursive(const double *llr, const int *info_mask,
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

        llr_left[i] = sign * min_val
                    + log1pexp_neg(a1 + a2)
                    - log1pexp_neg(fabs(a1 - a2));
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