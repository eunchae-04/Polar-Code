#include "plot.h"

#include <stdio.h>

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

void run_gnuplot_multiplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) return;

    fprintf(gp, "set terminal windows size 1300, 450\n");
    fprintf(gp, "set multiplot layout 1,3 title 'Polar Code Scenario Analyses (N=1024, K=512)' font ',13'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set logscale y 10\n");
    fprintf(gp, "set format y '10^{%%L}'\n");
    fprintf(gp, "unset yrange\n");
    fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");

    fprintf(gp, "set title '1. Baseline (Perfect)'\n");
    fprintf(gp,
            "plot 'simulation_baseline.txt' using 1:2 with linespoints lw 2 lc rgb 'purple' title 'BER', "
            "'simulation_baseline.txt' using 1:3 with linespoints lw 2 lc rgb 'cyan' title 'FER'\n");

    fprintf(gp, "set title '2. Scenario A (LLR Error, a=0.5)'\n");
    fprintf(gp,
            "plot 'simulation_scenario_a.txt' using 1:2 with linespoints lw 2 lc rgb 'red' title 'BER', "
            "'simulation_scenario_a.txt' using 1:3 with linespoints lw 2 lc rgb 'orange' title 'FER'\n");

    fprintf(gp, "set title '3. Scenario B (Design SNR=0dB)'\n");
    fprintf(gp,
            "plot 'simulation_scenario_b.txt' using 1:2 with linespoints lw 2 lc rgb 'blue' title 'BER', "
            "'simulation_scenario_b.txt' using 1:3 with linespoints lw 2 lc rgb 'dark-green' title 'FER'\n");

    fprintf(gp, "unset multiplot\n");
    PCLOSE(gp);
}