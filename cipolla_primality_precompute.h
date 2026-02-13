#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// interface to a multithreaded version for heavy computations
// -----------------------------------------------------------------------

#include "gmp.h"
#include <stdbool.h>
#include <stdint.h>

struct mod_precompute_t
{
    mpz_t a;           // Barrett coefficient 2^n32 mod m
    mpz_t b;           // Barrett coefficient 2^n32 div m
    mpz_t m;           // modulus
    mpz_t i;           // another integer
    uint64_t n;        // log2(modulus)
    uint64_t n2;       // Reduction threshold (n / 2) rounded down
    uint64_t n32;      // Barrett threshold (n * 3)/2, exactly n32 = n + n2
    bool special_case; // modulus has a special form
    bool montg;        // reduction requires montgomery form
    bool proth;        // modulus is proth number e * 2^n + 1
    bool power2me;     // modulus is 2^n - e
    bool power2pe;     // modulus is 2^n + e
    bool gmn;          // modulus is a * 2^n2 - b
    uint64_t e;        // small number part of the modulus
};

struct mod_precompute_t *mpz_mod_precompute(mpz_t n, bool verbose = false);
void mpz_mod_uncompute(mod_precompute_t *p);
void mpz_mod_fast_reduce(mpz_t r, mpz_t tmp, struct mod_precompute_t *p);
void mpz_mod_positive_reduce(mpz_t r, mpz_t tmp, struct mod_precompute_t *p);
void mpz_mod_div2(mpz_t r, struct mod_precompute_t *p);
void mpz_mod_to_montg(mpz_t v, struct mod_precompute_t *p);
void mpz_mod_from_montg(mpz_t v, mpz_t tmp, struct mod_precompute_t *p);
void mpz_mod_slow_reduce(mpz_t x, mpz_t m);
