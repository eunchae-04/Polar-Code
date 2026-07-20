#ifndef SIMULATION_H
#define SIMULATION_H

#include "polar_config.h"

void generate_info_bits(int *info_bits);
void build_u_from_mask(const int *mask, const int *info_bits, int *u);
void compute_llr(const double *rx, double sigma2_est, double *llr);
long count_bit_errors(const int *info_mask, const int *u_hat, const int *info_bits);
void run_simulation(const SimulationConfig *cfg, const int *fixed_mask);

#endif