/*
 * C 표준 rand() 함수 visual 결함 집중 검증 스크립트
 * 컴파일: gcc -O2 -o test_rand test_rand.c -lm
 * 실행: ./test_rand
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ZOOM_FILE "rand_zoom_grid.txt"
#define DATA_3D_FILE "rand_3d_planes.txt"

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

// Gnuplot Multiplot 시각화 (2D 미세 줌인 & 3D 평면 쏠림)
static void plot_visual_flaws(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) {
        printf("Notice: Gnuplot이 설치되어 있지 않거나 실행할 수 없습니다.\n");
        return;
    }

    fprintf(gp, "set terminal windows size 1200, 550\n");
    fprintf(gp, "set multiplot layout 1,2 title 'Visual Flaws of C Standard rand()' font ',13'\n");

    // [서브플롯 1] 2D 미세 구간 줌인 (양자화 바둑판 격자 확인)
    fprintf(gp, "set title '1. 2D Zoom-In [0.50 ~ 0.505] (Coarse Quantization Grid)'\n");
    fprintf(gp, "set xlabel 'u_i'\n");
    fprintf(gp, "set ylabel 'u_{i+1}'\n");
    fprintf(gp, "set xrange [0.500:0.505]\n");
    fprintf(gp, "set yrange [0.500:0.505]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '%s' using 1:2 with points pt 7 ps 0.4 lc rgb 'blue' title 'Discontinuous Grid'\n", ZOOM_FILE);

    // [서브플롯 2] 3D 입체 분포 (마사글리아 평면 쏠림 현상)
    fprintf(gp, "set title '2. 3D Space Pairs (Marsaglia Parallel Planes)'\n");
    fprintf(gp, "set xlabel 'u_i'\n");
    fprintf(gp, "set ylabel 'u_{i+1}'\n");
    fprintf(gp, "set zlabel 'u_{i+2}'\n");
    fprintf(gp, "set xrange [0:1]; set yrange [0:1]; set zrange [0:1]\n");
    fprintf(gp, "set view 55, 70\n"); // 평면 층이 잘 보이는 관측 각도 설정
    fprintf(gp, "set grid\n");
    fprintf(gp, "splot '%s' using 1:2:3 with points pt 7 ps 0.3 lc rgb 'red' title 'Parallel Sheets'\n", DATA_3D_FILE);

    fprintf(gp, "unset multiplot\n");
    PCLOSE(gp);
}

int main(void) {
    srand((unsigned)time(NULL));

    printf("=======================================================\n");
    printf("     C 표준 rand() 함수 시각적 결함 렌더링 중...        \n");
    printf("=======================================================\n");

    // 1. 미세 줌인용 영역 데이터 수집 ([0.500, 0.505] 상자 내부 점 수집)
    FILE *fp_zoom = fopen(ZOOM_FILE, "w");
    if (fp_zoom == NULL) return 1;

    int zoom_collected = 0;
    while (zoom_collected < 3000) {
        double u1 = (double)rand() / (double)RAND_MAX;
        double u2 = (double)rand() / (double)RAND_MAX;

        if (u1 >= 0.500 && u1 <= 0.505 && u2 >= 0.500 && u2 <= 0.505) {
            fprintf(fp_zoom, "%.7f %.7f\n", u1, u2);
            zoom_collected++;
        }
    }
    fclose(fp_zoom);

    // 2. 3D 입체 렌더링용 3연속 난수 튜플 수집 (u_i, u_{i+1}, u_{i+2})
    FILE *fp_3d = fopen(DATA_3D_FILE, "w");
    if (fp_3d == NULL) return 1;

    double u_i   = (double)rand() / (double)RAND_MAX;
    double u_i1  = (double)rand() / (double)RAND_MAX;

    for (int i = 0; i < 8000; i++) {
        double u_i2 = (double)rand() / (double)RAND_MAX;
        fprintf(fp_3d, "%.6f %.6f %.6f\n", u_i, u_i1, u_i2);
        u_i = u_i1;
        u_i1 = u_i2;
    }
    fclose(fp_3d);

    printf("데이터 수집 완료. Gnuplot 창을 띄웁니다.\n");
    printf("=======================================================\n");

    plot_visual_flaws();

    return 0;
}