#ifndef POLAR_CODEC_H
#define POLAR_CODEC_H

void polar_encode_recursive(const int *u, int *x, int n_len);
void polar_sc_decode_recursive(const double *llr, const int *info_mask,
                                int *u_hat, int *u_coded,
                                int n_len, int *mask_offset);

#endif