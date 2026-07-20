#include "polar_config.h"
#include "polar_math.h"
#include "simulation.h"
#include "plot.h"
#include "rng.h"

#include <stdio.h>

int main(void) {
    seed_xorshift64();

    SimulationConfig baseline = {
        .output_file = "simulation_baseline.txt",
        .title = "[1/3] Running Baseline Scenario...",
        .llr_scale = 1.0,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_a = {
        .output_file = "simulation_scenario_a.txt",
        .title = "[2/3] Running Scenario A (LLR Mismatch, alpha = 0.5)...",
        .llr_scale = 0.5,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_b = {
        .output_file = "simulation_scenario_b.txt",
        .title = "[3/3] Running Scenario B (Design SNR Mismatch = 0.0 dB Fixed)...",
        .llr_scale = 1.0,
        .use_fixed_mask = true,
        .fixed_mask_snr_db = 0.0
    };

    int fixed_info_mask[N];
    double sigma2_design = compute_sigma2_from_ebno_db(scenario_b.fixed_mask_snr_db);
    construct_frozen_mask(sigma2_design, fixed_info_mask);

    run_simulation(&baseline, NULL);
    run_simulation(&scenario_a, NULL);
    run_simulation(&scenario_b, fixed_info_mask);

    run_gnuplot_multiplot();

    printf("\nAll simulations completed! Gnuplot multiplot window popped up.\n");
    return 0;
}