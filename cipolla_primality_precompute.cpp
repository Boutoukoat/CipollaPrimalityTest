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

struct mod_precompute_t *mpz_mod_precompute(mpz_t n, bool verbose)
{
    mpz_t tmp;
    mod_precompute_t *p = (mod_precompute_t *)cipolla_allocate_function(sizeof(mod_precompute_t));

    p->special_case = false;
    p->proth = false;
    p->montg = false;
    p->power2pe = false;
    p->power2me = false;
    p->gmn = false;
    mpz_inits(tmp, p->a, p->b, p->m, p->inv, 0);
    mpz_inits(p->x_lo, p->x_hi, 0);
    p->n = mpz_sizeinbase(n, 2);
    p->n2 = 0;
    p->n32 = 0;
    p->e = 0;
    mpz_set(p->m, n);

    // check a power of 2 minus e
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
        // check a power of 2 plus e
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
        if (4 * --s > p->n * 3)
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
            if (mpz_cmp_ui(tmp, 1) == 0)
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
        // precompute a variant of Barrett reduction
        // b = 2^(3n/2) / n
        // a = 2^(3n/2) % n
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
        mpz_div_2exp(p->x_hi, r, p->n32 + p->n2);
        while (mpz_sgn(p->x_hi) != 0)
        {
            // (p->x_hi * a) << n/2 + p->x_lo
            mpz_mod_2exp(p->x_lo, r, p->n32 + p->n2);
            mpz_mul(tmp, p->x_hi, p->a);
            mpz_mul_2exp(p->x_hi, tmp, p->n2);
            mpz_add(r, p->x_hi, p->x_lo);
            mpz_div_2exp(p->x_hi, r, p->n32 + p->n2);
        }

        // reduce the number to approx 3*n/2 bits
        mpz_div_2exp(p->x_hi, r, p->n32);
        if (mpz_sgn(p->x_hi) != 0)
        {
            // p->x_hi * a + p->x_lo
            mpz_mod_2exp(p->x_lo, r, p->n32);
            mpz_mul(tmp, p->x_hi, p->a);
            mpz_add(r, tmp, p->x_lo);
        }

        // reduce the number to approx n bits (Barrett reduction)
        mpz_div_2exp(p->x_hi, r, p->n);
        mpz_mul(tmp, p->x_hi, p->b);
        mpz_div_2exp(p->x_hi, tmp, p->n2);
        mpz_mul(tmp, p->x_hi, p->m);
        mpz_sub(r, r, tmp);

        // reduce the number to exactly n bits
        mpz_div_2exp(tmp, r, p->n);
        while (mpz_sgn(tmp) != 0)
        {
            mpz_sub(r, r, p->m);
            mpz_div_2exp(tmp, r, p->n);
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

void mpz_mod_slow_reduce(mpz_t x, mpz_t m)
{
    // subtract the modulus until underflow
    while (mpz_cmp(x, m) >= 0)
    {
        mpz_sub(x, x, m);
    }
}
