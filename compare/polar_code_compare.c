/*
    gcc -O2 -Wall -o polar_code_compare polar_code_compare.c
    ./polar_code_compare
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(path) mkdir(path, 0755)
#endif

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

// 이미 만들어 둔 SC / SCL 결과 파일 경로 (compare 폴더 기준 상대경로)
#define SC_FILE  "../scd/result/simulation_results.txt"
#define SCL_FILE "../scl/result/simulation_results_scl.txt"

#define RESULT_DIR "result"
#define SUMMARY_FILE RESULT_DIR "/comparison_summary.txt"

#define MAX_POINTS 256

typedef struct {
    double snr_db;
    double ber;
    double fer;
} DataPoint;

static void ensure_result_dir(void) {
    if (MKDIR(RESULT_DIR) != 0 && errno != EEXIST) {
        printf("Warning: failed to create '%s' directory (errno=%d).\n", RESULT_DIR, errno);
    }
}

// "snr ber fer" 형식의 텍스트 파일을 읽어 배열에 담는다.
static int load_results(const char *path, DataPoint *points, int max_points) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        printf("Failed to open result file: %s\n", path);
        printf("  (compare 폴더 안에서 실행 중인지, scd/scl 폴더에 결과 파일이 있는지 확인해주세요.)\n");
        return -1;
    }

    int count = 0;
    while (count < max_points &&
           fscanf(fp, "%lf %lf %lf", &points[count].snr_db, &points[count].ber, &points[count].fer) == 3) {
        count++;
    }

    fclose(fp);
    return count;
}

// SC 대비 SCL이 SNR별로 몇 배 개선되었는지 콘솔에 출력하고 파일로도 저장한다.
static void print_comparison_summary(const DataPoint *sc, int sc_count,
                                      const DataPoint *scl, int scl_count) {
    int count = (sc_count < scl_count) ? sc_count : scl_count;

    FILE *fp = fopen(SUMMARY_FILE, "w");

    printf("\n========================================\n");
    printf(" SC vs SCL Comparison Summary\n");
    printf("========================================\n");
    printf("Eb/No(dB)  SC_BER        SCL_BER       BER_Gain(x)  SC_FER        SCL_FER       FER_Gain(x)\n");

    if (fp) {
        fprintf(fp, "Eb/No(dB) SC_BER SCL_BER BER_Gain SC_FER SCL_FER FER_Gain\n");
    }

    for (int i = 0; i < count; i++) {
        if (fabs(sc[i].snr_db - scl[i].snr_db) > 1e-6) {
            printf("Warning: SNR mismatch at index %d (SC=%.2f, SCL=%.2f). Skipping.\n",
                   i, sc[i].snr_db, scl[i].snr_db);
            continue;
        }

        double ber_gain = (scl[i].ber > 0.0) ? (sc[i].ber / scl[i].ber) : NAN;
        double fer_gain = (scl[i].fer > 0.0) ? (sc[i].fer / scl[i].fer) : NAN;

        printf("%6.2f     %.6e   %.6e   %9.2fx    %.6e   %.6e   %9.2fx\n",
               sc[i].snr_db, sc[i].ber, scl[i].ber, ber_gain,
               sc[i].fer, scl[i].fer, fer_gain);

        if (fp) {
            fprintf(fp, "%.2f %.6e %.6e %.4f %.6e %.6e %.4f\n",
                    sc[i].snr_db, sc[i].ber, scl[i].ber, ber_gain,
                    sc[i].fer, scl[i].fer, fer_gain);
        }
    }

    if (fp) {
        fclose(fp);
        printf("\nSaved to %s\n", SUMMARY_FILE);
    }
}

// gnuplot으로 SC/SCL의 BER, FER 4개 곡선을 한 그래프에 겹쳐 그린다.
static void run_gnuplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Warning: failed to launch gnuplot. Is it installed and in PATH?\n");
        return;
    }

    fprintf(gp, "set terminal windows size 900, 650\n");
    fprintf(gp, "set title 'Polar Code (N=1024, K=512): SC vs SCL(L=8) Performance'\n");
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'Eb/No (dB)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp,
            "plot '" SC_FILE  "' using 1:2 with linespoints lw 2 lc rgb 'red'    title 'SC BER', "
            "'" SC_FILE  "' using 1:3 with linespoints lw 2 lc rgb 'orange' title 'SC FER', "
            "'" SCL_FILE "' using 1:2 with linespoints lw 2 lc rgb 'blue'   title 'SCL BER', "
            "'" SCL_FILE "' using 1:3 with linespoints lw 2 lc rgb 'cyan'   title 'SCL FER'\n");

    PCLOSE(gp);
}

int main(void) {
    ensure_result_dir();

    DataPoint sc_points[MAX_POINTS];
    DataPoint scl_points[MAX_POINTS];

    int sc_count = load_results(SC_FILE, sc_points, MAX_POINTS);
    int scl_count = load_results(SCL_FILE, scl_points, MAX_POINTS);

    if (sc_count <= 0 || scl_count <= 0) {
        printf("결과 파일을 읽는 데 실패하여 종료합니다.\n");
        return 1;
    }

    printf("Loaded %d points from %s\n", sc_count, SC_FILE);
    printf("Loaded %d points from %s\n", scl_count, SCL_FILE);

    print_comparison_summary(sc_points, sc_count, scl_points, scl_count);

    run_gnuplot();
    printf("\nGnuplot window launched (if gnuplot is installed).\n");

    return 0;
}