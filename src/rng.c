#include "rng.h"
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t rng_state;

void seed_xorshift64(void) {
    rng_state = (uint64_t)time(NULL) ^ 0x5DEECE66DL;
    if (rng_state == 0) rng_state = 1;
}

uint64_t xorshift64(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

double rand_double(void) {
    return (double)xorshift64() / (double)UINT64_MAX;
}

double gauss_box_muller(void) {
    double u1, u2;
    do {
        u1 = rand_double();
    } while (u1 <= 1e-12);

    u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}