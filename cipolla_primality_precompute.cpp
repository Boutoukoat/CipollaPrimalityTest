#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "cipolla_primality_alloc.h"
#include "cipolla_primality_precompute.h"

typedef unsigned __int128 uint128_t;

// round up a to b boundary
static uint64_t uint64_round_up(uint64_t a, uint64_t b)
{
    a += b - 1;
    a /= b;
    a *= b;
    return a;
}

// Natural logarithm of a mpz_t number, with "double" accuracy +/- a few ulp
static double mpz_log(mpz_t n)
{
    uint64_t bl = mpz_sizeinbase(n, 2);
    double ll = 0.0;
    if (bl >= 64)
    {
        // compute the log from the 64 most significant bits, and scale by the amount shifted.
        mpz_t tmp;
        mpz_init(tmp);
        mpz_div_2exp(tmp, n, bl - 64);
        ll = log((double)mpz_get_ui(tmp)) + (bl - 64) * log((double)2.0);
        mpz_clear(tmp);
    }
    else
    {
        // compute the log of the only limb
        ll = log((double)mpz_get_ui(n));
    }
    return ll;
}

// check a number of (B^N-1)/(B-1) (sign = -1)
// check a number of (B^N+1)/(B+1) (sign = +1)
// return true if found
static bool mpz_isrepunit(mpz_t n, mpz_t tmp, uint64_t *pB, uint64_t *pN, int64_t *pSign)
{
    uint64_t B = 0;
    uint64_t N = 0;
    bool found = false;
    // try different values of B < 64 (arbitrarily)
    for (B = 3; B < 32; B++)
    {
        uint64_t B4 = mpz_mod_ui(tmp, n, B * B * B * B);
        if (B4 != ((B + 1) * B + 1) * B + 1)
        {
            // not a repunit signature
            continue;
        }
        // compute the log of B^N, i.e. find N
        mpz_mul_ui(tmp, n, B - 1);
        mpz_add_ui(tmp, tmp, 1);
        double ld = mpz_log(tmp);
        double lb = log((double)B);
        N = (uint64_t)floor(ld / lb + 0.5);
        if (N < 3)
        {
            // exponent is too small
            continue;
        }
        // rebuild the modulus to check the accuracy of B and N
        mpz_set_ui(tmp, B);
        mpz_pow_ui(tmp, tmp, N);
        mpz_sub_ui(tmp, tmp, 1);
        mpz_div_ui(tmp, tmp, B - 1);
        if (mpz_cmp(n, tmp) != 0)
        {
            // not the right number
            continue;
        }
        found = true;
        *pB = B;
        *pN = N;
        *pSign = -1;
        break;
    }
    return found;
}

struct mod_precompute_t *mpz_mod_precompute(mpz_t n, uint64_t log2a, bool verbose)
{
    mpz_t tmp;
    mod_precompute_t *p = (mod_precompute_t *)cipolla_allocate_function(sizeof(mod_precompute_t));

    p->special_case = false;
    p->proth = false;
    p->montg = false;
    p->power2pe = false;
    p->power2me = false;
    p->gmn = false;
    mpz_init(tmp);
    mpz_inits(p->a, p->b, p->m, p->inv, 0);
    mpz_inits(p->x_lo, p->x_hi, 0);

    p->n = mpz_sizeinbase(n, 2);
    p->n2 = 0;
    p->n32 = 0;
    p->e = 0;
    mpz_set(p->m, n);

    // check a power of 2 minus e     (e < 2^64)
    mpz_set_ui(tmp, 1);
    mpz_mul_2exp(tmp, tmp, p->n);
    mpz_sub(tmp, tmp, n); // tmp = 2^n - modulus
    p->power2me = (p->n > 128 && mpz_sgn(tmp) >= 0 && mpz_size(tmp) <= 1);
    if (p->power2me)
    {
        p->e = mpz_get_ui(tmp);
        p->special_case = true;
    }

    if (!p->special_case)
    {
        // check a power of 2 plus e    (e < 2^64)
        mpz_set_ui(tmp, 1);
        mpz_mul_2exp(tmp, tmp, p->n - 1);
        mpz_sub(tmp, n, tmp); // tmp = modulus - 2^n
        p->power2pe = (p->n > 128 && mpz_sgn(tmp) >= 0 && mpz_size(tmp) <= 1);
        if (p->power2pe)
        {
            p->e = mpz_get_ui(tmp);
            p->special_case = true;
        }
    }

    if (!p->special_case)
    {
        // check a Proth number  b * 2^n2 + 1
        p->n2 = (p->n + 1) / 2;
        mpz_mod_2exp(tmp, n, p->n2);
        if (p->n2 >= 64 && mpz_cmp_ui(tmp, 1) == 0)
        {
            while (mpz_tstbit(n, p->n2) == 0)
            {
                p->n2 += 1;
            }
            // the least significant half of the number is 0x0001....1
            // precompute the most significant part of the number
            mpz_div_2exp(p->b, n, p->n2);
            p->n32 = p->n + (p->n >> 1);
            mpz_set_ui(tmp, 1);
            mpz_mul_2exp(tmp, tmp, p->n32);
            mpz_mod(p->a, tmp, n);
            p->proth = true;
            p->montg = true;
            p->special_case = true;
        }
    }

    if (!p->special_case)
    {
        // check a generalized mersenne number a * 2^n2 - b
        // start from the middle of the modulus
        uint64_t s = p->n / 2;
        while (mpz_tstbit(n, s) == 1)
        {
            s += 1;
        }
        // make sure reduction is worth it  (TODO)
	if (3 * --s > p-> n)
	{
            mpz_set_ui(tmp, 1);
            mpz_mul_2exp(tmp, tmp, s);
            mpz_mod_2exp(p->b, n, s);
            mpz_sub(p->b, tmp, p->b);
            mpz_div_2exp(p->a, n, s);
            mpz_add_ui(p->a, p->a, 1);
            while ((mpz_get_ui(p->a) & 1) == 0)
            {
                mpz_div_2exp(p->a, p->a, 1);
                s += 1;
            }
            mpz_gcd(tmp, p->a, n);
            if (mpz_cmp_ui(tmp, 1) == 0 && 2 * (mpz_sizeinbase(p->a, 2) + mpz_sizeinbase(p->b, 2)) < p->n)
            {
                mpz_mul(tmp, p->a, p->a);
                mpz_invert(p->inv, tmp, n);
                p->n2 = s;
                p->montg = true;
                p->gmn = true;
                p->special_case = true;
            }
	}
    }

    if (!p->special_case)
    {
        // detect a repunit (B^N - 1)/(B-1)
        // detect a repunit (B^N + 1)/(B+1)     N odd
        //
        // (repunit is prime iff N is prime)
        //
        // e.g. (10^15-1)/(10-1) = 111111111111111
        // e.g. (10^15+1)/(10+1) = 90909090909091
        //
        uint64_t B = 0;
        uint64_t N = 0;
        int64_t sign = 0;
        if (mpz_isrepunit(p->m, tmp, &B, &N, &sign))
        {
            // p->repunit = true;
            // p->special_case = true;
            if (verbose)
            {
                printf("Repunit (%lu^%lu%s1)/(%lu%s1)\n", B, N, (sign > 0 ? "+" : "-"), B, (sign > 0 ? "+" : "-"));
                printf("Not yet optimized\n");
                // TODO : ifigure out a modular reduction mod B^N
            }
        }
    }

    if (!p->special_case)
    {
        // precompute a variant of Barrett reduction mod n where numbers to reduce are O(6 * a *n * n)
        // b = 2^(3n/2) / n
        // a = 2^(3n/2) % n
        p->n = uint64_round_up(p->n + log2a + 3, 2); // fine tuning of the number of bits
        p->n2 = p->n >> 1;
        p->n32 = p->n + p->n2;
        mpz_set_ui(tmp, 1);
        mpz_mul_2exp(tmp, tmp, p->n32);
        mpz_divmod(p->b, p->a, tmp, n);
    }

    mpz_clears(tmp, 0);
    if (verbose)
    {
        if (p->power2pe)
        {
            printf("Modular reduction optimized for numbers 2^s + e\n");
        }
        if (p->power2me)
        {
            printf("Modular reduction optimized for numbers 2^s - e\n");
        }
        if (p->proth)
        {
            printf("Modular reduction optimized for numbers e*2^s + 1\n");
        }
        if (p->gmn)
        {
            printf("Modular reduction optimized for numbers a*2^s - b\n");
        }
        if (!p->special_case)
        {
            printf("Modular reduction not optimized\n");
        }
    }
    return p;
}

void mpz_mod_uncompute(mod_precompute_t *p)
{
    if (p)
    {
        mpz_clears(p->a, p->b, p->m, p->inv, 0);
        mpz_clears(p->x_lo, p->x_hi, 0);
        memset(p, 0, sizeof(struct mod_precompute_t));
        cipolla_free_function(p, sizeof(struct mod_precompute_t));
    }
}

// input
//     r   : a number to reduce, can be much larger than modulus^2 by magnitude orders
//     tmp : scratch area
//     p   : precomputed constants and flags where p->m is the modulus
// output
//     r   : a reduced number >= 0 and < 2*modulus
void mpz_mod_fast_reduce(mpz_t r, mpz_t tmp, struct mod_precompute_t *p)
{
    if (p->special_case)
    {
        // special reduction for modulus = b * 2^n + 1
        if (p->proth)
        {
            if (mpz_sizeinbase(r, 2) > 2 * p->n + 2)
            {
                mpz_div_2exp(p->x_hi, r, p->n32);
                if (mpz_sgn(p->x_hi) != 0)
                {
                    // p->x_hi * a + p->x_lo
                    mpz_mod_2exp(p->x_lo, r, p->n32);
                    mpz_mul(tmp, p->x_hi, p->a);
                    mpz_add(r, tmp, p->x_lo);
                }
            }

            // reduce the number to approx n bits (Montgomery reduction)
            mpz_div_2exp(p->x_hi, r, p->n2);
            mpz_mod_2exp(p->x_lo, r, p->n2);
            mpz_mul(tmp, p->x_lo, p->b);
            mpz_sub(tmp, tmp, p->x_hi);
            mpz_div_2exp(p->x_hi, tmp, p->n2);
            mpz_mod_2exp(p->x_lo, tmp, p->n2);
            mpz_mul(tmp, p->x_lo, p->b);
            mpz_sub(r, tmp, p->x_hi);
            if (mpz_sgn(r) < 0)
            {
                mpz_add(r, r, p->m);
            }
            else if (mpz_cmp(r, p->m) >= 0)
            {
                mpz_sub(r, r, p->m);
            }
        }

        // special reduction for modulus = 2^n - e
        else if (p->power2me)
        {
            // while (hi != 0) r = lo + hi * e
            mpz_div_2exp(p->x_hi, r, p->n);
            while (mpz_sgn(p->x_hi) != 0)
            {
                mpz_mod_2exp(p->x_lo, r, p->n);
                mpz_mul_ui(tmp, p->x_hi, p->e);
                mpz_add(r, p->x_lo, tmp);
                mpz_div_2exp(p->x_hi, r, p->n);
            }
        }

        // special reduction for modulus = 2^n + e
        else if (p->power2pe)
        {
            // while (hi != 0) r = lo - hi * e
            mpz_div_2exp(p->x_hi, r, p->n - 1);
            while (mpz_cmp_ui(p->x_hi, 1) > 0)
            {
                mpz_mod_2exp(p->x_lo, r, p->n - 1);
                mpz_mul_ui(tmp, p->x_hi, p->e);
                if (mpz_cmp(p->x_lo, tmp) >= 0)
                {
                    //    lo - hi * e
                    mpz_sub(r, p->x_lo, tmp);
                }
                else
                {
                    //    (lo + k * m) - hi * e
                    mpz_set(p->x_hi, tmp);
                    mpz_div_2exp(tmp, p->x_hi, p->n - 1);
                    mpz_add_ui(tmp, tmp, 1);
                    mpz_mul(r, tmp, p->m);
                    mpz_add(r, r, p->x_lo);
                    mpz_sub(r, r, p->x_hi);
                }
                mpz_div_2exp(p->x_hi, r, p->n - 1);
            }
        }
        else if (p->gmn)
        {
            // special reduction for modulus = a*2^n2 - b for a, b small
            mpz_div_2exp(p->x_hi, r, p->n2);
            mpz_mod_2exp(p->x_lo, r, p->n2);
            mpz_mul(p->x_hi, p->x_hi, p->b);
            mpz_mul(p->x_lo, p->x_lo, p->a);
            mpz_add(r, p->x_lo, p->x_hi);
            mpz_div_2exp(p->x_hi, r, p->n2);
            mpz_mod_2exp(p->x_lo, r, p->n2);
            mpz_mul(p->x_hi, p->x_hi, p->b);
            mpz_mul(p->x_lo, p->x_lo, p->a);
            mpz_add(r, p->x_lo, p->x_hi);
        }
    }
    else
    {
        // reduce the number to approx 2*n bits
        // mpz_div_2exp(p->x_hi, r, p->n32 + p->n2);
        // while (mpz_sgn(p->x_hi) != 0)
        while (mpz_sizeinbase(r, 2) > p->n32 + p->n2 + 2)
        {
            // (p->x_hi * a) << n/2 + p->x_lo
            mpz_div_2exp(p->x_hi, r, p->n32 + p->n2);
            mpz_mod_2exp(p->x_lo, r, p->n32 + p->n2);
            mpz_mul(tmp, p->x_hi, p->a);
            mpz_mul_2exp(p->x_hi, tmp, p->n2);
            mpz_add(r, p->x_hi, p->x_lo);
            //    mpz_div_2exp(p->x_hi, r, p->n32 + p->n2);
        }

        // reduce the number to approx 3*n/2 bits
        if (mpz_sizeinbase(r, 2) > p->n32 + 2)
        {
            // p->x_hi * a + p->x_lo
            mpz_div_2exp(p->x_hi, r, p->n32);
            mpz_mod_2exp(p->x_lo, r, p->n32);
            mpz_mul(tmp, p->x_hi, p->a);
            mpz_add(r, p->x_lo, tmp);
        }

        // reduce the number to approx n bits (Barrett reduction)
        mpz_div_2exp(p->x_hi, r, p->n);
        mpz_mul(tmp, p->x_hi, p->b);
        mpz_div_2exp(p->x_hi, tmp, p->n2);
        mpz_mul(tmp, p->x_hi, p->m);
        mpz_sub(r, r, tmp);

        // reduce the number to exactly n bits
        while (mpz_sizeinbase(r, 2) > p->n + 2)
        {
            mpz_sub(r, r, p->m);
        }
    }
}

void mpz_mod_positive_reduce(mpz_t r, mpz_t tmp, struct mod_precompute_t *p)
{
    // check r < 0
    if (mpz_sgn(r) < 0)
    {
        // make number positive by adding a large multiple of the modulus
        int bits = mpz_sizeinbase(r, 2) - p->n;
        bits = bits < 0 ? 1 : bits + 1;
        mpz_mul_2exp(tmp, p->m, bits);
        mpz_add(r, r, tmp);
    }
}

void mpz_mod_div2(mpz_t r, struct mod_precompute_t *p)
{
    // check r odd
    if (mpz_get_ui(r) & 1)
    {
        // make number even by adding the modulus
        mpz_add(r, r, p->m);
    }
    // divide an even number by 2
    mpz_div_2exp(r, r, 1);
}

void mpz_mod_to_montg(mpz_t v, struct mod_precompute_t *p)
{
    if (p->montg)
    {
        if (p->proth)
        {
            mpz_mul_2exp(v, v, 2 * p->n2);
            mpz_mod(v, v, p->m);
        }
        else if (p->gmn)
        {
            mpz_mul(v, v, p->inv);
            mpz_mod(v, v, p->m);
        }
    }
}

void mpz_mod_from_montg(mpz_t v, mpz_t tmp, struct mod_precompute_t *p)
{
    if (p->montg)
    {
        mpz_mod_fast_reduce(v, tmp, p);
    }
}

void mpz_mod_slow_reduce(mpz_t x, struct mod_precompute_t *p)
{
    if (!p->special_case)
    {
        mpz_mod(x, x, p->m);
    }
    else
    {
        // subtract the modulus until underflow
        while (mpz_cmp(x, p->m) >= 0)
        {
            mpz_sub(x, x, p->m);
        }
    }
}
