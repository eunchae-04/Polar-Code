#ifndef RNG_H
#define RNG_H

#include <stdint.h>

void seed_xorshift64(void);
uint64_t xorshift64(void);
double rand_double(void);
double gauss_box_muller(void);

#endif