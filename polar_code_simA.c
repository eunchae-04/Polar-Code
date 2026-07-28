#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 1024
#define K 512
#define R ((double)K / (double)N)

#define EBNO_START 0.0
#define EBNO_END   3.0
#define EBNO_STEP  0.25
#define TARGET_ERRORS 10000

// (EBNO_END - EBNO_START) / EBNO_STEP + 1 = (3.0 - 0.0) / 0.25 + 1 = 13
// EBNO_START/END/STEP을 바꾸면 이 값도 함께 맞춰줘야 합니다.
#define NUM_SNR_POINTS 13

#define DEFAULT_LENA_INPUT "image.png"
#define DEFAULT_IMAGE_SNR_DB 3.0
#define RESULT_DIR "result"

// Scenario A alpha sweep 설정: 실제 채널을 고정하고 alpha(=LLR 스케일 오차 계수)를
// 여러 값으로 훑으며 BER/이미지 결과를 관찰. 0<alpha<1(LLR 확대)과 alpha>1(LLR 축소)을
// 모두 포함하도록 대표값을 선정.
#define ALPHA_SWEEP_TRUE_SNR_DB DEFAULT_IMAGE_SNR_DB
static const double ALPHA_SWEEP_VALUES[] = {0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0};
#define ALPHA_SWEEP_COUNT ((int)(sizeof(ALPHA_SWEEP_VALUES) / sizeof(ALPHA_SWEEP_VALUES[0])))

// Scenario B design-SNR sweep 설정: 실제 채널을 이 SNR로 고정하고,
// 마스크 설계 SNR을 DESIGN_SNR_SWEEP_START~END까지 훑으며 BER 변화를 관찰.
#define DESIGN_SNR_SWEEP_TRUE_DB 3.0
#define DESIGN_SNR_SWEEP_START (-2.0)
#define DESIGN_SNR_SWEEP_END   5.0
#define DESIGN_SNR_SWEEP_STEP  0.5

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(path) mkdir(path, 0755)
#endif
#include <errno.h>

// result/ 폴더가 없으면 생성. 이미 있으면(EEXIST) 조용히 통과.
static void ensure_result_dir(void) {
    if (MKDIR(RESULT_DIR) != 0 && errno != EEXIST) {
        printf("Warning: failed to create '%s' directory (errno=%d). Output files may fail to save.\n",
               RESULT_DIR, errno);
    }
}

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#define PYTHON_CMD "C:/msys64/ucrt64/bin/python.exe"
#else
#define POPEN popen
#define PCLOSE pclose
#define PYTHON_CMD "python3"
#endif

#define TMP_INPUT_PGM  "__polar_tmp_input.pgm"
#define TMP_OUTPUT_PGM "__polar_tmp_output.pgm"

typedef struct {
    const char *output_file;
    const char *title;
    double llr_scale;             // sigma2_est = sigma2_true * llr_scale
    bool use_fixed_mask;          // true면 fixed_mask_snr_db로 생성한 마스크를 고정 사용
    double fixed_mask_snr_db;     // fixed mask 생성용 SNR
} SimulationConfig;

typedef struct {
    int width;
    int height;
    unsigned char *pixels;
} GrayImage;

typedef struct {
    long pixel_differences;
    double mean_absolute_error;
    double mean_squared_error;
    double psnr;
    unsigned char min_difference;
    unsigned char max_difference;
} ImageComparisonStats;

// =================================================================
// 난수 생성기
// =================================================================
static uint64_t rng_state;

static void seed_xorshift64(void) {
    rng_state = (uint64_t)time(NULL) ^ 0x5DEECE66DL;
    if (rng_state == 0) rng_state = 1;
}

static uint64_t xorshift64(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static double rand_double(void) {
    return (double)xorshift64() / (double)UINT64_MAX;
}

static double gauss_box_muller(void) {
    double u1, u2;
    do {
        u1 = rand_double();
    } while (u1 <= 1e-12);
    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// =================================================================
// Polar Code 구성
// =================================================================
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

static double compute_sigma2_from_ebno_db(double ebno_db) {
    double ebno_linear = pow(10.0, ebno_db / 10.0);
    return 1.0 / (2.0 * R * ebno_linear);
}

static void construct_frozen_mask(double sigma2, int *info_mask) {
    double mu[N];
    mu[0] = 2.0 / sigma2;

    int n = 10; // N = 2^10 = 1024
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

// =================================================================
// Encode / Decode
// =================================================================
static void polar_encode_recursive(const int *u, int *x, int n_len) {
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

static void polar_sc_decode_recursive(const double *llr, const int *info_mask,
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
                    + log(1.0 + exp(-(a1 + a2)))
                    - log(1.0 + exp(-fabs(a1 - a2)));
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

// =================================================================
// 공통 유틸리티
// =================================================================
static void generate_info_bits(int *info_bits) {
    for (int i = 0; i < K; i++) {
        info_bits[i] = (int)(xorshift64() & 1ULL);
    }
}

static void build_u_from_mask(const int *mask, const int *info_bits, int *u) {
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

static void compute_llr(const double *rx, double sigma2_est, double *llr) {
    for (int i = 0; i < N; i++) {
        llr[i] = (2.0 / sigma2_est) * rx[i];
    }
}

static long count_bit_errors(const int *info_mask, const int *u_hat, const int *info_bits) {
    int info_ptr = 0;
    long bit_errors = 0;

    for (int i = 0; i < N; i++) {
        if (info_mask[i]) {
            if (u_hat[i] != info_bits[info_ptr++]) bit_errors++;
        }
    }
    return bit_errors;
}

static void free_gray_image(GrayImage *image) {
    if (image->pixels != NULL) {
        free(image->pixels);
        image->pixels = NULL;
    }
    image->width = 0;
    image->height = 0;
}

static bool read_pgm_token(FILE *fp, char *buffer, size_t buffer_size) {
    int ch;

    do {
        ch = fgetc(fp);
        if (ch == '#') {
            while (ch != '\n' && ch != EOF) ch = fgetc(fp);
        }
    } while (ch != EOF && isspace((unsigned char)ch));

    if (ch == EOF) return false;

    size_t len = 0;
    while (ch != EOF && !isspace((unsigned char)ch) && ch != '#') {
        if (len + 1 < buffer_size) {
            buffer[len++] = (char)ch;
        }
        ch = fgetc(fp);
    }

    if (ch == '#') {
        while (ch != '\n' && ch != EOF) ch = fgetc(fp);
    }

    buffer[len] = '\0';
    return len > 0;
}

static bool read_pgm_image(const char *file_path, GrayImage *image) {
    FILE *fp = fopen(file_path, "rb");
    if (fp == NULL) {
        printf("Failed to open input image: %s\n", file_path);
        return false;
    }

    int max_val = 0;
    bool ok = false;

    if (fscanf(fp, " P5 %d %d %d", &image->width, &image->height, &max_val) != 3) {
        printf("Failed to parse PGM header: %s\n", file_path);
        goto cleanup;
    }
    if (image->width <= 0 || image->height <= 0 || max_val != 255) {
        printf("Unsupported PGM dimensions or max value: %d x %d, max=%d\n",
               image->width, image->height, max_val);
        goto cleanup;
    }

    int separator = fgetc(fp);
    if (separator == EOF) {
        printf("PGM header ended unexpectedly: %s\n", file_path);
        goto cleanup;
    }

    size_t pixel_count = (size_t)image->width * (size_t)image->height;
    image->pixels = (unsigned char *)malloc(pixel_count);
    if (image->pixels == NULL) {
        printf("Failed to allocate %zu bytes for image pixels.\n", pixel_count);
        goto cleanup;
    }

    if (fread(image->pixels, 1, pixel_count, fp) != pixel_count) {
        printf("Failed to read %zu image pixels from %s\n", pixel_count, file_path);
        goto cleanup;
    }

    ok = true;

cleanup:
    if (!ok) {
        free_gray_image(image);
    }
    fclose(fp);
    return ok;
}

static bool has_extension_ci(const char *path, const char *extension) {
    size_t path_len = strlen(path);
    size_t ext_len = strlen(extension);
    if (path_len < ext_len) return false;

    const char *path_ext = path + (path_len - ext_len);
    for (size_t i = 0; i < ext_len; i++) {
        unsigned char a = (unsigned char)path_ext[i];
        unsigned char b = (unsigned char)extension[i];
        if (tolower(a) != tolower(b)) return false;
    }
    return true;
}

static bool convert_image_to_pgm(const char *input_path, const char *output_pgm_path) {
    char command[2048];
    int written = snprintf(
        command,
        sizeof(command),
        "%s -c \"from PIL import Image; import sys; Image.open(sys.argv[1]).convert('L').save(sys.argv[2])\" \"%s\" \"%s\"",
        PYTHON_CMD,
        input_path,
        output_pgm_path
    );

    if (written < 0 || (size_t)written >= sizeof(command)) {
        printf("Failed to build image conversion command.\n");
        return false;
    }

    return system(command) == 0;
}

static void build_preview_path(const char *output_path, char *preview_path, size_t preview_size) {
    const char *dot = strrchr(output_path, '.');
    if (dot != NULL && has_extension_ci(dot, ".png")) {
        size_t base_len = (size_t)(dot - output_path);
        snprintf(preview_path, preview_size, "%.*s_preview.png", (int)base_len, output_path);
        return;
    }

    snprintf(preview_path, preview_size, "%s_preview.png", output_path);
}

static bool create_image_preview(const char *input_path, const char *decoded_path, const char *preview_path) {
    char command[4096];
    int written = snprintf(
        command,
        sizeof(command),
        "%s -c \"from PIL import Image, ImageChops, ImageOps, ImageDraw; import sys; src=Image.open(sys.argv[1]).convert('L'); dec=Image.open(sys.argv[2]).convert('L'); diff=ImageOps.autocontrast(ImageChops.difference(src, dec)); pad=12; label_h=28; w=max(src.width, dec.width, diff.width); h=max(src.height, dec.height, diff.height); canvas=Image.new('L', (w*3 + pad*4, h + label_h + pad*2), 255); d=ImageDraw.Draw(canvas); canvas.paste(src, (pad, pad + label_h)); d.text((pad, pad), 'Original', fill=0); canvas.paste(dec, (w + pad*2, pad + label_h)); d.text((w + pad*2, pad), 'Decoded', fill=0); canvas.paste(diff, (w*2 + pad*3, pad + label_h)); d.text((w*2 + pad*3, pad), 'Diff', fill=0); canvas.save(sys.argv[3])\" \"%s\" \"%s\" \"%s\"",
        PYTHON_CMD,
        input_path,
        decoded_path,
        preview_path
    );

    if (written < 0 || (size_t)written >= sizeof(command)) {
        printf("Failed to build preview image command.\n");
        return false;
    }

    return system(command) == 0;
}

static void open_image_file(const char *path) {
#ifdef _WIN32
    char command[2048];
    if (snprintf(command, sizeof(command), "cmd /c start \"\" \"%s\"", path) > 0) {
        system(command);
    }
#else
    (void)path;
#endif
}

static bool load_input_image_any_format(const char *input_path, GrayImage *image, char *temp_pgm_path, size_t temp_pgm_size) {
    if (has_extension_ci(input_path, ".pgm")) {
        return read_pgm_image(input_path, image);
    }

    if (snprintf(temp_pgm_path, temp_pgm_size, "%s", TMP_INPUT_PGM) < 0) {
        return false;
    }

    if (!convert_image_to_pgm(input_path, temp_pgm_path)) {
        printf("Failed to convert input image to PGM: %s\n", input_path);
        return false;
    }

    if (!read_pgm_image(temp_pgm_path, image)) {
        remove(temp_pgm_path);
        return false;
    }

    return true;
}

static bool convert_pgm_to_requested_output(const char *pgm_path, const char *output_path) {
    if (has_extension_ci(output_path, ".pgm")) {
        remove(output_path);
        return rename(pgm_path, output_path) == 0;
    }

    if (!convert_image_to_pgm(pgm_path, output_path)) {
        printf("Failed to convert output PGM to target format: %s\n", output_path);
        return false;
    }

    return true;
}

static bool write_pgm_image(const char *file_path, const GrayImage *image) {
    FILE *fp = fopen(file_path, "wb");
    if (fp == NULL) return false;

    size_t pixel_count = (size_t)image->width * (size_t)image->height;
    bool ok = fprintf(fp, "P5\n%d %d\n255\n", image->width, image->height) > 0
           && fwrite(image->pixels, 1, pixel_count, fp) == pixel_count;

    fclose(fp);
    return ok;
}

static void bytes_to_bits_msb(const unsigned char *bytes, size_t byte_count, int *bits) {
    for (size_t i = 0; i < byte_count; i++) {
        unsigned char value = bytes[i];
        for (int bit = 0; bit < 8; bit++) {
            bits[i * 8 + bit] = (value >> (7 - bit)) & 1;
        }
    }
}

static void bits_to_bytes_msb(const int *bits, size_t bit_count, unsigned char *bytes) {
    size_t byte_count = bit_count / 8;
    for (size_t i = 0; i < byte_count; i++) {
        unsigned char value = 0;
        for (int bit = 0; bit < 8; bit++) {
            value = (unsigned char)((value << 1) | (bits[i * 8 + bit] & 1));
        }
        bytes[i] = value;
    }
}

static ImageComparisonStats compare_gray_images(const GrayImage *before, const GrayImage *after) {
    ImageComparisonStats stats = {0, 0.0, 0.0, 0.0, 255, 0};
    size_t pixel_count = (size_t)before->width * (size_t)before->height;

    for (size_t i = 0; i < pixel_count; i++) {
        int diff = (int)before->pixels[i] - (int)after->pixels[i];
        int abs_diff = abs(diff);

        if (abs_diff != 0) {
            stats.pixel_differences++;
        }

        stats.mean_absolute_error += (double)abs_diff;
        stats.mean_squared_error += (double)(diff * diff);

        if ((unsigned char)abs_diff < stats.min_difference) {
            stats.min_difference = (unsigned char)abs_diff;
        }
        if ((unsigned char)abs_diff > stats.max_difference) {
            stats.max_difference = (unsigned char)abs_diff;
        }
    }

    if (pixel_count > 0) {
        stats.mean_absolute_error /= (double)pixel_count;
        stats.mean_squared_error /= (double)pixel_count;
    }

    if (stats.mean_squared_error > 0.0) {
        stats.psnr = 10.0 * log10((255.0 * 255.0) / stats.mean_squared_error);
    } else {
        stats.psnr = INFINITY;
    }

    return stats;
}

static void transmit_image_bits(const int *input_bits, size_t input_bit_count,
                                double snr_db, bool use_fixed_mask,
                                double fixed_mask_snr_db, double llr_scale,
                                int *output_bits) {
    double sigma2_true = compute_sigma2_from_ebno_db(snr_db);
    double sigma2_est = sigma2_true * llr_scale;

    int fixed_mask[N];
    int *info_mask_source = NULL;
    if (use_fixed_mask) {
        double sigma2_design = compute_sigma2_from_ebno_db(fixed_mask_snr_db);
        construct_frozen_mask(sigma2_design, fixed_mask);
        info_mask_source = fixed_mask;
    }

    size_t block_count = (input_bit_count + K - 1) / K;
    size_t padded_bit_count = block_count * (size_t)K;

    for (size_t block = 0; block < block_count; block++) {
        int info_bits[K] = {0};
        int u[N];
        int x[N];
        double rx[N];
        double llr[N];
        int u_hat[N];
        int u_coded[N];
        int info_mask[N];
        int mask_offset = 0;

        if (use_fixed_mask) {
            for (int i = 0; i < N; i++) info_mask[i] = info_mask_source[i];
        } else {
            construct_frozen_mask(sigma2_true, info_mask);
        }

        size_t block_start = block * (size_t)K;
        size_t remaining = (input_bit_count > block_start) ? (input_bit_count - block_start) : 0;
        size_t copy_count = (remaining < (size_t)K) ? remaining : (size_t)K;

        for (size_t i = 0; i < copy_count; i++) {
            info_bits[i] = input_bits[block_start + i];
        }

        build_u_from_mask(info_mask, info_bits, u);
        polar_encode_recursive(u, x, N);
        awgn_channel(x, sigma2_true, rx);
        compute_llr(rx, sigma2_est, llr);
        polar_sc_decode_recursive(llr, info_mask, u_hat, u_coded, N, &mask_offset);

        size_t output_offset = block * (size_t)K;
        size_t out_ptr = 0;
        for (int i = 0; i < N; i++) {
            if (info_mask[i]) {
                output_bits[output_offset + out_ptr] = u_hat[i];
                out_ptr++;
            }
        }

        (void)padded_bit_count;
    }
}

static void run_image_mode(const char *scenario_label, const char *input_path, const char *output_path,
                            double snr_db, bool use_fixed_mask, double fixed_mask_snr_db, double llr_scale) {
    GrayImage input_image = {0, 0, NULL};
    GrayImage output_image = {0, 0, NULL};
    char temp_input_pgm[64] = {0};
    char temp_output_pgm[64] = TMP_OUTPUT_PGM;
    char preview_path[512] = {0};

    printf("\n--- [%s] Image Transmission ---\n", scenario_label);

    if (!load_input_image_any_format(input_path, &input_image, temp_input_pgm, sizeof(temp_input_pgm))) {
        printf("Failed to read input image: %s\n", input_path);
        return;
    }

    size_t pixel_count = (size_t)input_image.width * (size_t)input_image.height;
    size_t input_bit_count = pixel_count * 8;
    size_t block_count = (input_bit_count + K - 1) / K;
    size_t output_bit_count = block_count * (size_t)K;

    int *input_bits = (int *)calloc(output_bit_count, sizeof(int));
    int *output_bits = (int *)calloc(output_bit_count, sizeof(int));
    if (input_bits == NULL || output_bits == NULL) {
        printf("Failed to allocate bit buffers for image transmission.\n");
        free(input_bits);
        free(output_bits);
        free_gray_image(&input_image);
        return;
    }

    bytes_to_bits_msb(input_image.pixels, pixel_count, input_bits);
    transmit_image_bits(input_bits, input_bit_count, snr_db, use_fixed_mask, fixed_mask_snr_db, llr_scale, output_bits);

    output_image.width = input_image.width;
    output_image.height = input_image.height;
    output_image.pixels = (unsigned char *)malloc(pixel_count);
    if (output_image.pixels == NULL) {
        printf("Failed to allocate output image buffer.\n");
        free(input_bits);
        free(output_bits);
        free_gray_image(&input_image);
        return;
    }

    bits_to_bytes_msb(output_bits, input_bit_count, output_image.pixels);

    ImageComparisonStats stats = compare_gray_images(&input_image, &output_image);

    long bit_errors = 0;
    for (size_t i = 0; i < input_bit_count; i++) {
        if (input_bits[i] != output_bits[i]) bit_errors++;
    }

    if (!write_pgm_image(temp_output_pgm, &output_image)) {
        printf("Failed to write output image: %s\n", output_path);
    } else {
        bool output_ready = convert_pgm_to_requested_output(temp_output_pgm, output_path);
        remove(temp_output_pgm);

        if (!output_ready) {
            printf("Failed to finalize output image: %s\n", output_path);
        } else {
            build_preview_path(output_path, preview_path, sizeof(preview_path));
            if (create_image_preview(input_path, output_path, preview_path)) {
                printf("Preview image: %s\n", preview_path);
                open_image_file(preview_path);
            }

            printf("Scenario: %s\n", scenario_label);
            printf("Input: %s\n", input_path);
            printf("Output: %s\n", output_path);
            printf("SNR: %.2f dB\n", snr_db);
            printf("Bit errors: %ld / %zu\n", bit_errors, input_bit_count);
            printf("Pixel differences: %ld / %zu\n", stats.pixel_differences, pixel_count);
            printf("MAE: %.6f\n", stats.mean_absolute_error);
            printf("MSE: %.6f\n", stats.mean_squared_error);
            if (isfinite(stats.psnr)) {
                printf("PSNR: %.6f dB\n", stats.psnr);
            } else {
                printf("PSNR: infinite (no pixel changes)\n");
            }
            printf("Diff range: %u to %u\n", stats.min_difference, stats.max_difference);
        }
    }

    if (temp_input_pgm[0] != '\0') {
        remove(temp_input_pgm);
    }
    free(input_bits);
    free(output_bits);
    free_gray_image(&input_image);
    free_gray_image(&output_image);
}

static void run_simulation(const SimulationConfig *cfg, const int *fixed_mask,
                            double *snr_out, double *ber_out, double *fer_out) {
    FILE *fp = fopen(cfg->output_file, "w");
    if (fp == NULL) {
        printf("Failed to open output file: %s\n", cfg->output_file);
        return;
    }

    printf("\n%s\n", cfg->title);
    printf("Eb/No (dB)\tBER\t\tFER\n");

    int point_idx = 0;
    for (double snr_db = EBNO_START; snr_db <= EBNO_END + 1e-9; snr_db += EBNO_STEP) {
        double sigma2_true = compute_sigma2_from_ebno_db(snr_db);
        double sigma2_est = sigma2_true * cfg->llr_scale;

        int info_mask[N];
        if (cfg->use_fixed_mask) {
            for (int i = 0; i < N; i++) info_mask[i] = fixed_mask[i];
        } else {
            construct_frozen_mask(sigma2_est, info_mask);
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

        if (point_idx < NUM_SNR_POINTS) {
            if (snr_out) snr_out[point_idx] = snr_db;
            if (ber_out) ber_out[point_idx] = final_ber;
            if (fer_out) fer_out[point_idx] = final_fer;
        }
        point_idx++;
    }

    fclose(fp);
}

// baseline 대비 scenario A/B의 BER/FER이 SNR별로 몇 배 나빠지는지 요약 출력 + 파일 저장.
// 교수님이 물어본 "알파 오추정이 얼마나 차이 나느냐"에 숫자로 답하기 위한 함수.
static void print_scenario_comparison(const double *snr, int count,
                                       const double *base_ber, const double *base_fer,
                                       const double *a_ber, const double *a_fer,
                                       const double *b_ber, const double *b_fer) {
    FILE *fp = fopen(RESULT_DIR "/scenario_comparison_summary.txt", "w");

    printf("\n========================================\n");
    printf(" Scenario Comparison Summary (vs Baseline)\n");
    printf("========================================\n");
    printf("Eb/No(dB)  Base_BER      ScenA_BER     A/Base(x)   ScenB_BER     B/Base(x)\n");

    if (fp) {
        fprintf(fp, "Eb/No(dB) Base_BER ScenA_BER A_over_Base ScenB_BER B_over_Base\n");
    }

    for (int i = 0; i < count; i++) {
        double ratio_a = (base_ber[i] > 0.0) ? (a_ber[i] / base_ber[i]) : NAN;
        double ratio_b = (base_ber[i] > 0.0) ? (b_ber[i] / base_ber[i]) : NAN;

        printf("%6.2f     %.6e   %.6e   %9.2fx   %.6e   %9.2fx\n",
               snr[i], base_ber[i], a_ber[i], ratio_a, b_ber[i], ratio_b);

        if (fp) {
            fprintf(fp, "%.2f %.6e %.6e %.4f %.6e %.4f\n",
                    snr[i], base_ber[i], a_ber[i], ratio_a, b_ber[i], ratio_b);
        }

        (void)base_fer; (void)a_fer; (void)b_fer;
    }

    if (fp) {
        printf("\nSaved to %s/scenario_comparison_summary.txt\n", RESULT_DIR);
        fclose(fp);
    }
}

static void run_gnuplot_alpha_sweep(const char *data_file, double true_snr) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) return;

    fprintf(gp, "set terminal windows size 900, 500\n");
    fprintf(gp, "set title 'Scenario A: BER/FER vs alpha (true SNR = %.2f dB)'\n", true_snr);
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale x\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'alpha (sigma2_est / sigma2_true)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp, "set arrow from 1, graph 0 to 1, graph 1 nohead lc rgb 'gray' dt 2\n");
    fprintf(gp, "set label 'alpha=1 (baseline)' at 1, graph 0.95 offset 1,0 tc rgb 'gray'\n");
    fprintf(gp,
            "plot '%s' using 1:2 with linespoints lw 2 pt 7 lc rgb 'red' title 'BER', "
            "'%s' using 1:3 with linespoints lw 2 pt 7 lc rgb 'blue' title 'FER'\n",
            data_file, data_file);

    PCLOSE(gp);
}

static void run_gnuplot_design_snr_sweep(const char *data_file, double true_snr) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (gp == NULL) return;

    fprintf(gp, "set terminal windows size 900, 500\n");
    fprintf(gp, "set title 'Scenario B: BER/FER vs Design SNR (true channel = %.2f dB)'\n", true_snr);
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'Design SNR (dB)'\n");
    fprintf(gp, "set ylabel 'Error Probability'\n");
    fprintf(gp, "set arrow from %.2f, graph 0 to %.2f, graph 1 nohead lc rgb 'gray' dt 2\n", true_snr, true_snr);
    fprintf(gp, "set label 'design=true (ideal)' at %.2f, graph 0.95 offset 1,0 tc rgb 'gray'\n", true_snr);
    fprintf(gp,
            "plot '%s' using 1:2 with linespoints lw 2 pt 7 lc rgb 'red' title 'BER', "
            "'%s' using 1:3 with linespoints lw 2 pt 7 lc rgb 'blue' title 'FER'\n",
            data_file, data_file);

    PCLOSE(gp);
}
static void run_gnuplot_multiplot(void) {
    FILE *gp = POPEN("gnuplot -persist", "w");
    fprintf(gp, "set multiplot layout 1,3 title 'Polar Code Scenario Analyses (N=1024, K=512)' font ',13'\n");
    fprintf(gp, "set datafile separator whitespace\n");
    fprintf(gp, "set logscale y\n");
    fprintf(gp, "set grid\n");
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

static void print_usage(const char *program_name) {
    printf("Usage:\n");
    printf("  %s                       # run the 3 built-in BER/FER simulations + per-scenario image tests\n", program_name);
    printf("  %s --alpha <value>       # ad hoc alpha (sigma2_est = alpha * sigma2_true) test vs baseline\n", program_name);
    printf("  %s --sweep-a [true_snr]  # Scenario A alpha sweep (BER + image per alpha) at a fixed true SNR (default %.2f dB)\n",
           program_name, ALPHA_SWEEP_TRUE_SNR_DB);
    printf("  %s --sweep-b [true_snr]  # Scenario B design-SNR sweep at a fixed true channel SNR (default %.2f dB)\n",
           program_name, DESIGN_SNR_SWEEP_TRUE_DB);
    printf("  %s --image input.pgm output.pgm [snr_db]\n", program_name);
    printf("\n");
    printf("Image mode supports binary grayscale PGM (P5) and, via Pillow, common formats (png/bmp/jpg).\n");
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--sweep-a") == 0) {
        double true_snr = ALPHA_SWEEP_TRUE_SNR_DB;
        if (argc >= 3) {
            true_snr = atof(argv[2]);
        }

        seed_xorshift64();
        ensure_result_dir();

        double sigma2_true = compute_sigma2_from_ebno_db(true_snr);

        // Scenario A 정의: 마스크는 항상 실제 채널 기준으로 정확하게 구성 (LLR 스케일만 오차)
        int correct_mask[N];
        construct_frozen_mask(sigma2_true, correct_mask);

        double alpha_ber[ALPHA_SWEEP_COUNT];
        double alpha_fer[ALPHA_SWEEP_COUNT];

        for (int ai = 0; ai < ALPHA_SWEEP_COUNT; ai++) {
            double alpha = ALPHA_SWEEP_VALUES[ai];
            double sigma2_est = sigma2_true * alpha;

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
                build_u_from_mask(correct_mask, info_bits, u);
                polar_encode_recursive(u, x, N);
                awgn_channel(x, sigma2_true, rx);
                compute_llr(rx, sigma2_est, llr);
                polar_sc_decode_recursive(llr, correct_mask, u_hat, u_coded, N, &mask_offset);

                long bit_errors = count_bit_errors(correct_mask, u_hat, info_bits);
                if (bit_errors > 0) total_frame_errors++;
                total_errors += bit_errors;
                total_bits += K;
                total_frames++;
            }

            alpha_ber[ai] = (double)total_errors / (double)total_bits;
            alpha_fer[ai] = (double)total_frame_errors / (double)total_frames;

            printf("  [%d/%d] alpha = %.2f done (BER = %.6e)\n", ai + 1, ALPHA_SWEEP_COUNT, alpha, alpha_ber[ai]);

            char img_out[256];
            snprintf(img_out, sizeof(img_out), RESULT_DIR "/lena_output_alpha_sweep_%.2f.png", alpha);
            char img_label[64];
            snprintf(img_label, sizeof(img_label), "Alpha Sweep alpha=%.2f", alpha);

            run_image_mode(img_label, DEFAULT_LENA_INPUT, img_out, true_snr, false, 0.0, alpha);
        }

        // alpha=1.0 지점을 baseline 기준으로 찾아 배율(ratio) 계산
        double baseline_ber = -1.0;
        for (int ai = 0; ai < ALPHA_SWEEP_COUNT; ai++) {
            if (fabs(ALPHA_SWEEP_VALUES[ai] - 1.0) < 1e-9) {
                baseline_ber = alpha_ber[ai];
                break;
            }
        }

        printf("\n========================================\n");
        printf(" Scenario A Alpha Sensitivity Sweep (true channel = %.2f dB)\n", true_snr);
        printf("========================================\n");
        printf("Alpha    BER            FER            Ratio_to_alpha1   Note\n");

        char sweep_file[256];
        snprintf(sweep_file, sizeof(sweep_file), RESULT_DIR "/scenario_a_alpha_sweep_true%.2f.txt", true_snr);
        FILE *fp = fopen(sweep_file, "w");
        if (fp) fprintf(fp, "Alpha BER FER Ratio_to_alpha1\n");

        for (int ai = 0; ai < ALPHA_SWEEP_COUNT; ai++) {
            double alpha = ALPHA_SWEEP_VALUES[ai];
            double ratio = (baseline_ber > 0.0) ? (alpha_ber[ai] / baseline_ber) : NAN;
            const char *note = (fabs(alpha - 1.0) < 1e-9) ? "<- baseline (alpha=1)" : "";

            printf("%5.2f    %.6e   %.6e   %13.2fx     %s\n",
                   alpha, alpha_ber[ai], alpha_fer[ai], ratio, note);
            if (fp) fprintf(fp, "%.2f %.6e %.6e %.4f\n", alpha, alpha_ber[ai], alpha_fer[ai], ratio);
        }

        if (fp) {
            fclose(fp);
            printf("\nSaved to %s\n", sweep_file);
        }

        run_gnuplot_alpha_sweep(sweep_file, true_snr);

        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "--sweep-b") == 0) {
        double true_snr = DESIGN_SNR_SWEEP_TRUE_DB;
        if (argc >= 3) {
            true_snr = atof(argv[2]);
        }

        seed_xorshift64();
        ensure_result_dir();

        double sigma2_true = compute_sigma2_from_ebno_db(true_snr);

        int num_points = (int)((DESIGN_SNR_SWEEP_END - DESIGN_SNR_SWEEP_START) / DESIGN_SNR_SWEEP_STEP + 1.5);
        double *design_snr_arr = (double *)malloc(sizeof(double) * (size_t)num_points);
        double *design_ber_arr = (double *)malloc(sizeof(double) * (size_t)num_points);
        double *design_fer_arr = (double *)malloc(sizeof(double) * (size_t)num_points);
        int count = 0;

        for (double design_snr = DESIGN_SNR_SWEEP_START; design_snr <= DESIGN_SNR_SWEEP_END + 1e-9; design_snr += DESIGN_SNR_SWEEP_STEP) {
            int fixed_mask[N];
            double sigma2_design = compute_sigma2_from_ebno_db(design_snr);
            construct_frozen_mask(sigma2_design, fixed_mask);

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
                build_u_from_mask(fixed_mask, info_bits, u);
                polar_encode_recursive(u, x, N);
                awgn_channel(x, sigma2_true, rx);
                compute_llr(rx, sigma2_true, llr); // LLR 스케일은 정확 (마스크만 미스매치)
                polar_sc_decode_recursive(llr, fixed_mask, u_hat, u_coded, N, &mask_offset);

                long bit_errors = count_bit_errors(fixed_mask, u_hat, info_bits);
                if (bit_errors > 0) total_frame_errors++;
                total_errors += bit_errors;
                total_bits += K;
                total_frames++;
            }

            if (count < num_points) {
                design_snr_arr[count] = design_snr;
                design_ber_arr[count] = (double)total_errors / (double)total_bits;
                design_fer_arr[count] = (double)total_frame_errors / (double)total_frames;
            }
            printf("  [%d/%d] design_snr = %.2f dB done (BER = %.6e)\n",
                   count + 1, num_points, design_snr, design_ber_arr[count]);
            count++;
        }

        double matched_ber = -1.0;
        for (int i = 0; i < count; i++) {
            if (fabs(design_snr_arr[i] - true_snr) < DESIGN_SNR_SWEEP_STEP / 2.0) {
                matched_ber = design_ber_arr[i];
                break;
            }
        }

        printf("\n========================================\n");
        printf(" Scenario B Design-SNR Sensitivity Sweep (true channel = %.2f dB)\n", true_snr);
        printf("========================================\n");
        printf("Design_SNR(dB)   BER            FER            Ratio_to_matched   Note\n");

        char sweep_file[256];
        snprintf(sweep_file, sizeof(sweep_file), RESULT_DIR "/scenario_b_design_snr_sweep_true%.2f.txt", true_snr);
        FILE *fp = fopen(sweep_file, "w");
        if (fp) fprintf(fp, "Design_SNR(dB) BER FER Ratio_to_matched\n");

        for (int i = 0; i < count; i++) {
            double ratio = (matched_ber > 0.0) ? (design_ber_arr[i] / matched_ber) : NAN;
            const char *note = (fabs(design_snr_arr[i] - true_snr) < DESIGN_SNR_SWEEP_STEP / 2.0) ? "<- matched (ideal)" : "";

            printf("%8.2f       %.6e   %.6e   %15.2fx   %s\n",
                   design_snr_arr[i], design_ber_arr[i], design_fer_arr[i], ratio, note);
            if (fp) fprintf(fp, "%.2f %.6e %.6e %.4f\n",
                             design_snr_arr[i], design_ber_arr[i], design_fer_arr[i], ratio);
        }

        if (fp) {
            fclose(fp);
            printf("\nSaved to %s\n", sweep_file);
        }

        free(design_snr_arr);
        free(design_ber_arr);
        free(design_fer_arr);

        run_gnuplot_design_snr_sweep(sweep_file, true_snr);

        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "--alpha") == 0) {
        if (argc < 3) {
            printf("Usage: %s --alpha <value>   (e.g. --alpha 2.0 for alpha>1, --alpha 0.5 for alpha<1)\n", argv[0]);
            return 1;
        }

        double alpha = atof(argv[2]);
        if (alpha <= 0.0) {
            printf("Alpha must be a positive number.\n");
            return 1;
        }

        seed_xorshift64();
        ensure_result_dir();

        char ber_file[256];
        snprintf(ber_file, sizeof(ber_file), RESULT_DIR "/simulation_custom_alpha_%.2f.txt", alpha);
        char baseline_file[256];
        snprintf(baseline_file, sizeof(baseline_file), RESULT_DIR "/simulation_baseline_for_alpha_%.2f.txt", alpha);
        char title_buf[128];
        snprintf(title_buf, sizeof(title_buf), "[2/2] Custom Scenario A (alpha = %.2f)...", alpha);

        SimulationConfig custom_baseline = {
            .output_file = baseline_file,
            .title = "[1/2] Baseline (for alpha comparison)...",
            .llr_scale = 1.0,
            .use_fixed_mask = false,
            .fixed_mask_snr_db = 0.0
        };

        SimulationConfig custom_alpha = {
            .output_file = ber_file,
            .title = title_buf,
            .llr_scale = alpha,
            .use_fixed_mask = false,
            .fixed_mask_snr_db = 0.0
        };

        double snr_points[NUM_SNR_POINTS];
        double base_ber[NUM_SNR_POINTS], base_fer[NUM_SNR_POINTS];
        double custom_ber[NUM_SNR_POINTS], custom_fer[NUM_SNR_POINTS];

        run_simulation(&custom_baseline, NULL, snr_points, base_ber, base_fer);
        run_simulation(&custom_alpha, NULL, NULL, custom_ber, custom_fer);

        char summary_file[256];
        snprintf(summary_file, sizeof(summary_file), RESULT_DIR "/comparison_alpha_%.2f.txt", alpha);
        FILE *sfp = fopen(summary_file, "w");

        printf("\n========================================\n");
        printf(" alpha = %.2f vs Baseline\n", alpha);
        printf("========================================\n");
        printf("Eb/No(dB)  Base_BER      Custom_BER    Ratio(x)\n");
        if (sfp) fprintf(sfp, "Eb/No(dB) Base_BER Custom_BER Ratio\n");

        for (int i = 0; i < NUM_SNR_POINTS; i++) {
            double ratio = (base_ber[i] > 0.0) ? (custom_ber[i] / base_ber[i]) : NAN;
            printf("%6.2f     %.6e   %.6e   %8.2fx\n", snr_points[i], base_ber[i], custom_ber[i], ratio);
            if (sfp) fprintf(sfp, "%.2f %.6e %.6e %.4f\n", snr_points[i], base_ber[i], custom_ber[i], ratio);
        }

        if (sfp) {
            fclose(sfp);
            printf("\nSaved to %s\n", summary_file);
        }

        char img_out[256];
        snprintf(img_out, sizeof(img_out), RESULT_DIR "/lena_output_alpha_%.2f.png", alpha);
        char img_label[64];
        snprintf(img_label, sizeof(img_label), "Custom alpha=%.2f", alpha);

        run_image_mode(img_label, DEFAULT_LENA_INPUT, img_out,
                       DEFAULT_IMAGE_SNR_DB, false, 0.0, alpha);

        return 0;
    }

    if (argc >= 2 && strcmp(argv[1], "--image") == 0) {
        if (argc < 4) {
            print_usage(argv[0]);
            return 1;
        }

        double snr_db = 3.0;
        if (argc >= 5) {
            snr_db = atof(argv[4]);
        }

        seed_xorshift64();
        run_image_mode("Custom", argv[2], argv[3], snr_db, false, 0.0, 1.0);
        return 0;
    }

    seed_xorshift64();
    ensure_result_dir();

    SimulationConfig baseline = {
        .output_file = RESULT_DIR "/simulation_baseline.txt",
        .title = "[1/3] Running Baseline Scenario...",
        .llr_scale = 1.0,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_a = {
        .output_file = RESULT_DIR "/simulation_scenario_a.txt",
        .title = "[2/3] Running Scenario A (LLR Mismatch, alpha = 0.5)...",
        .llr_scale = 0.5,
        .use_fixed_mask = false,
        .fixed_mask_snr_db = 0.0
    };

    SimulationConfig scenario_b = {
        .output_file = RESULT_DIR "/simulation_scenario_b.txt",
        .title = "[3/3] Running Scenario B (Design SNR Mismatch = 0.0 dB Fixed)...",
        .llr_scale = 1.0,
        .use_fixed_mask = true,
        .fixed_mask_snr_db = 0.0
    };

    int fixed_info_mask[N];
    double sigma2_design = compute_sigma2_from_ebno_db(scenario_b.fixed_mask_snr_db);
    construct_frozen_mask(sigma2_design, fixed_info_mask);

    double snr_points[NUM_SNR_POINTS];
    double baseline_ber[NUM_SNR_POINTS], baseline_fer[NUM_SNR_POINTS];
    double scenario_a_ber[NUM_SNR_POINTS], scenario_a_fer[NUM_SNR_POINTS];
    double scenario_b_ber[NUM_SNR_POINTS], scenario_b_fer[NUM_SNR_POINTS];

    run_simulation(&baseline, NULL, snr_points, baseline_ber, baseline_fer);
    run_simulation(&scenario_a, NULL, NULL, scenario_a_ber, scenario_a_fer);
    run_simulation(&scenario_b, fixed_info_mask, NULL, scenario_b_ber, scenario_b_fer);

    print_scenario_comparison(snr_points, NUM_SNR_POINTS,
                               baseline_ber, baseline_fer,
                               scenario_a_ber, scenario_a_fer,
                               scenario_b_ber, scenario_b_fer);

    printf("\n========================================\n");
    printf(" Image Transmission Test per Scenario (%s)\n", DEFAULT_LENA_INPUT);
    printf("========================================\n");

    run_image_mode("Baseline", DEFAULT_LENA_INPUT, RESULT_DIR "/lena_output_baseline.png",
                   DEFAULT_IMAGE_SNR_DB, baseline.use_fixed_mask, baseline.fixed_mask_snr_db, baseline.llr_scale);

    run_image_mode("Scenario A (LLR Mismatch)", DEFAULT_LENA_INPUT, RESULT_DIR "/lena_output_scenario_a.png",
                   DEFAULT_IMAGE_SNR_DB, scenario_a.use_fixed_mask, scenario_a.fixed_mask_snr_db, scenario_a.llr_scale);

    run_image_mode("Scenario B (Design SNR Mismatch)", DEFAULT_LENA_INPUT, RESULT_DIR "/lena_output_scenario_b.png",
                   DEFAULT_IMAGE_SNR_DB, scenario_b.use_fixed_mask, scenario_b.fixed_mask_snr_db, scenario_b.llr_scale);

    run_gnuplot_multiplot();

    printf("\nAll simulations completed! Gnuplot multiplot window popped up.\n");
    return 0;
}