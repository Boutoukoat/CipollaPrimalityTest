// -----------------------------------------------------------------------
// Primality check based on Cipolla algorithm
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Pari/gp code :
//
// isCipollaPrime(n)={
// if((n==2||n==3||n==5||n==7||n==11),return(1));
// if(n<13||n%2==0,return(0));
// s=sqrtint(n)+1;d=(s*s)%n;
// while(ispower(d,2),s=s+1;d=(s*s)%n;);
// t=s+1;t2=(t*t-d)%n;j=kronecker(t2,n);if(j==0,return(0));
// while(j!=-1,t=t+1;t2=(t*t-d)%n;j=kronecker(t2,n);if(j==0,return(0)));
// A=Mod(Mod(x+t,n),x^2-t2);B=A^((n+1)/2);S=B*B;
// if(S==d,return(1),return(0));
// }
//
// for(i=100,10000000,n=2*i+1;if(isCipollaPrime(n) != isprime(n),print(n)))
//
// -----------------------------------------------------------------------

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "cipolla_primality.h"
#include "cipolla_primality_alloc.h"
#include "cipolla_primality_precompute.h"

typedef unsigned __int128 uint128_t;

// x % (2^b -1)
static uint64_t mpz_mod_mersenne(mpz_t x, uint64_t b)
{
    uint64_t mask = (1ull << b) - 1;
    unsigned s = mpz_size(x);
    const mp_limb_t *array = mpz_limbs_read(x);
    uint128_t t, r = 0;

    if ((b & (b - 1)) == 0)
    {
        // b is an exact power of 2
        for (unsigned i = 0; i < s; i += 1)
        {
            r += array[i];
        }
    }
    else
    {
        // add, shift first, and reduce after
        unsigned c = 0;
        for (unsigned i = 0; i < b; i += 1)
        {
            t = 0;
            for (unsigned j = i; j < s; j += b)
            {
                t += array[j];
            }
            t = (t & mask) + (t >> b);
            t = t << c;
            r += t;
            r = (r & mask) + (r >> b);
            c += 64;
            while (c >= b)
            {
                c -= b;
            }
        }
    }
    // reduce mod 2^b - 1 to a 64 bit number
    while (r >> 64)
    {
        r = (r & mask) + (r >> b);
    }
    return (uint64_t)r;
}

// (u << 64 + v) mod n
static inline uint64_t longmod(uint64_t u, uint64_t v, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r, a;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(u), "1"(v), "r"(n));
    return r;
#else
    uint128_t t = ((uint128_t)u << 64) + v;
    return t % n;
#endif
}

// (u) mod n
static inline uint64_t longlongmod(uint128_t u, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r = (uint64_t)(u >> 64), a = (uint64_t)u;
    asm("divq %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    return u % n;
#endif
}

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %3" : "=d"(r), "=a"(a) : "1"(a), "r"(b));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t tmp = (uint128_t)a * b;
    tmp %= n;
    return tmp;
#endif
}

static inline uint64_t squaremod(uint64_t a, uint64_t n)
{
#ifdef __x86_64__
    uint64_t r;
    asm("mul %2" : "=d"(r), "=a"(a) : "1"(a));
    asm("div %4" : "=d"(r), "=a"(a) : "0"(r), "1"(a), "r"(n));
    return r;
#else
    uint128_t tmp = (uint128_t)a * a;
    tmp %= n;
    return tmp;
#endif
}

// (u << s) mod n
static inline uint64_t shiftmod(uint64_t u, uint64_t s, uint64_t n)
{
#ifdef __x86_164__
    uint64_t r;
    asm("xor %0, %0\n shldq %b3, %1, %0\n shlxq %3, %1, %%rax\n divq %2"
        : "=&d"(r)
        : "r"(u), "r"(n), "c"(s)
        : "flags", "%rax");
    return r;
#else
    uint128_t t = (uint128_t)u;
    t <<= s;
    return t % n;
#endif
}

// count leading zeroed bits
static inline uint64_t uint64_lzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("lzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_clzll(a);
#endif
}

// count trailing zeroed bits
static inline uint64_t uint64_tzcnt(uint64_t a)
{
#ifdef __x86_64__
    uint64_t r;
    asm("tzcnt %1,%0" : "=r"(r) : "r"(a));
    return r;
#else
    return __builtin_ctzll(a);
#endif
}

// cipolla floor(log_2) function for integers
// log(1) = 0
// log(2) = 1
// log(3) = 1
// log(4) = 2 ...
static inline uint64_t uint64_log_2(uint64_t a)
{
    return 63 - uint64_lzcnt(a);
}

typedef enum sieve_e
{
    COMPOSITE_FOR_SURE,
    PRIME_FOR_SURE,
    UNDECIDED
} sieve_t;

// sieve small primes
static sieve_t uint64_composite_sieve(uint64_t a)
{
    if (a <= 152)
    {
        // return COMPOSITE_FOR_SURE for small numbers which are composite for sure , without checking further.
        sieve_t stooopid_prime_table[] = {
            PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE,
            COMPOSITE_FOR_SURE, COMPOSITE_FOR_SURE, PRIME_FOR_SURE,     COMPOSITE_FOR_SURE, PRIME_FOR_SURE,
            COMPOSITE_FOR_SURE};
        return stooopid_prime_table[a];
    }

    // divisibility is based on Barrett modular reductions by constants, use modular multiplications
    if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
        return COMPOSITE_FOR_SURE; // divisible by 3
    if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
        return COMPOSITE_FOR_SURE; // divisible by 5
    if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
        return COMPOSITE_FOR_SURE; // divisible by 7
    if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
        return COMPOSITE_FOR_SURE; // divisible by 11
    if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
        return COMPOSITE_FOR_SURE; // divisible by 13
    if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
        return COMPOSITE_FOR_SURE; // divisible by 17
    if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
        return COMPOSITE_FOR_SURE; // divisible by 19
    if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
        return COMPOSITE_FOR_SURE; // divisible by 23
    if ((uint64_t)(a * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
        return COMPOSITE_FOR_SURE; // divisible by 29
    if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
        return COMPOSITE_FOR_SURE; // divisible by 31
    if (a < 37 * 37)
        return PRIME_FOR_SURE;
    if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
        return COMPOSITE_FOR_SURE; // divisible by 37
    if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
        return COMPOSITE_FOR_SURE; // divisible by 41
    if ((uint64_t)(a * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
        return COMPOSITE_FOR_SURE; // divisible by 43
    if ((uint64_t)(a * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
        return COMPOSITE_FOR_SURE; // divisible by 47
    if ((uint64_t)(a * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
        return COMPOSITE_FOR_SURE; // divisible by 53
    if ((uint64_t)(a * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
        return COMPOSITE_FOR_SURE; // divisible by 59
    if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
        return COMPOSITE_FOR_SURE; // divisible by 61
    if ((uint64_t)(a * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
        return COMPOSITE_FOR_SURE; // divisible by 67
    if ((uint64_t)(a * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
        return COMPOSITE_FOR_SURE; // divisible by 71
    if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
        return COMPOSITE_FOR_SURE; // divisible by 73
    if ((uint64_t)(a * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
        return COMPOSITE_FOR_SURE; // divisible by 79
    if ((uint64_t)(a * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
        return COMPOSITE_FOR_SURE; // divisible by 83
    if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
        return COMPOSITE_FOR_SURE; // divisible by 89
    if ((uint64_t)(a * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
        return COMPOSITE_FOR_SURE; // divisible by 97
    if (a < 101 * 101)
        return PRIME_FOR_SURE;
    if ((uint64_t)(a * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
        return COMPOSITE_FOR_SURE; // divisible by 101
    if ((uint64_t)(a * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
        return COMPOSITE_FOR_SURE; // divisible by 103
    if ((uint64_t)(a * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
        return COMPOSITE_FOR_SURE; // divisible by 107
    if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
        return COMPOSITE_FOR_SURE; // divisible by 109
    if ((uint64_t)(a * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
        return COMPOSITE_FOR_SURE; // divisible by 113
    if ((uint64_t)(a * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
        return COMPOSITE_FOR_SURE; // divisible by 127
    if ((uint64_t)(a * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
        return COMPOSITE_FOR_SURE; // divisible by 131
    if ((uint64_t)(a * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
        return COMPOSITE_FOR_SURE; // divisible by 137
    if ((uint64_t)(a * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
        return COMPOSITE_FOR_SURE; // divisible by 139
    if ((uint64_t)(a * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
        return COMPOSITE_FOR_SURE; // divisible by 149
    if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
        return COMPOSITE_FOR_SURE; // divisible by 151

    // no factor less than 157
    if (a < 157 * 157)
        return PRIME_FOR_SURE; // prime
    return UNDECIDED;
}

static int uint64_jacobi(uint64_t x, uint64_t y)
{
    // assert((y & 1) == 1);
    if (y == 1 || x == 1)
    {
        return 1;
    }

    if (x == 2)
    {
        // char j[4] = { -1,-1,1,1};
        // return j[((y - 3) >> 1) % 4];
        return ((y + 2) & 4) ? -1 : 1;
    }
    if (x == 3)
    {
        char j[6] = {0, (char)-1, (char)-1, 0, 1, 1};
        return j[((y - 3) >> 1) % 6];
    }
    if (x == 5)
    {
        char j[5] = {(char)-1, 0, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 5];
    }
    if (x == 7)
    {
        char j[14] = {1, (char)-1, 0, 1, (char)-1, (char)-1, (char)-1, (char)-1, 1, 0, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 14];
    }
    if (x == 11)
    {
        char j[22] = {(char)-1, 1,        1,        1, 0, (char)-1, (char)-1, (char)-1, 1, (char)-1, (char)-1, 1,
                      (char)-1, (char)-1, (char)-1, 0, 1, 1,        1,        (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 22];
    }
    if (x == 13)
    {
        char j[13] = {1, (char)-1, (char)-1, 1, (char)-1, 0, (char)-1, 1, (char)-1, (char)-1, 1, 1, 1};
        return j[((y - 3) >> 1) % 13];
    }
    if (x == 17)
    {
        char j[17] = {(char)-1, (char)-1, (char)-1, 1,        (char)-1, 1,        1, 0, 1,
                      1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1, 1};
        return j[((y - 3) >> 1) % 17];
    }
    if (x == 19)
    {
        char j[19] = {1,        1, (char)-1, 1,        (char)-1, (char)-1, 1,        1,        0,       (char)-1,
                      (char)-1, 1, 1,        (char)-1, 1,        (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 38;
        return t >= 19 ? -j[t - 19] : j[t];
    }
    if (x == 23)
    {
        char j[23] = {(char)-1, (char)-1, 1,        1, 1,        1,        1,        (char)-1,
                      1,        (char)-1, 0,        1, (char)-1, 1,        (char)-1, (char)-1,
                      (char)-1, (char)-1, (char)-1, 1, 1,        (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 46;
        return t >= 23 ? -j[t - 23] : j[t];
    }
    if (x == 29)
    {
        char j[29] = {(char)-1, 1, 1,        1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1,
                      1,        1, (char)-1, 0, (char)-1, 1, 1,        (char)-1, (char)-1, (char)-1,
                      (char)-1, 1, (char)-1, 1, 1,        1, (char)-1, 1,        1};
        return j[((y - 3) >> 1) % 29];
    }
    if (x == 31)
    {
        char j[31] = {1,        1, (char)-1, 1,        1, (char)-1, 1,        (char)-1, (char)-1, (char)-1, 1,
                      1,        1, (char)-1, 0,        1, (char)-1, (char)-1, (char)-1, 1,        1,        1,
                      (char)-1, 1, (char)-1, (char)-1, 1, (char)-1, (char)-1, (char)-1, (char)-1};
        unsigned t = ((y - 3) >> 1) % 62;
        return t >= 31 ? -j[t - 31] : j[t];
    }

    int t = 1;
    uint64_t a = x;
    uint64_t n = y;
    unsigned v = n & 7;
    unsigned c = (v == 3) || (v == 5);
    while (a)
    {
        v = __builtin_ctzll(a);
        a >>= v;
        t = (c & (v & 1)) ? -t : t;

        if (a < n)
        {
            uint64_t tmp = a;
            a = n;
            n = tmp;
            t = ((a & n & 3) == 3) ? -t : t;
            v = n & 7;
            c = (v == 3) || (v == 5);
        }

        a -= n;
    }

    return (n == 1) ? t : 0;
}

static sieve_t mpz_composite_sieve(mpz_t n)
{
    // detect small 64-bit numbers
    unsigned sz = mpz_sizeinbase(n, 2);
    if (sz <= 64)
    {
        uint64_t a = mpz_get_ui(n);
        sieve_t sv = uint64_composite_sieve(a);
        return sv;
    }
    else
    {
        // large number, do the modular reduction in 2 steps
        // step 1 :
        //   reduce by a multiple of small factors
        // step 2:
        //    divisibility is based on Barrett modular reductions by constants, use modular multiplications

        uint64_t a;

        // 2^60-1 is divisible by 3,5,7,11,13,31,41,61,151 ...
        a = mpz_mod_mersenne(n, 60);
        if ((uint64_t)(a * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
            return COMPOSITE_FOR_SURE; // divisible by 3
        if ((uint64_t)(a * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
            return COMPOSITE_FOR_SURE; // divisible by 5
        if ((uint64_t)(a * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
            return COMPOSITE_FOR_SURE; // divisible by 7
        if ((uint64_t)(a * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
            return COMPOSITE_FOR_SURE; // divisible by 11
        if ((uint64_t)(a * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
            return COMPOSITE_FOR_SURE; // divisible by 13
        if ((uint64_t)(a * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
            return COMPOSITE_FOR_SURE; // divisible by 31
        if ((uint64_t)(a * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
            return COMPOSITE_FOR_SURE; // divisible by 41
        if ((uint64_t)(a * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
            return COMPOSITE_FOR_SURE; // divisible by 61
        if ((uint64_t)(a * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
            return COMPOSITE_FOR_SURE; // divisible by 151

        // 2^56-1 is divisible by 3, 5, 17, 29, 43, 113, 127, .....
        a = mpz_mod_mersenne(n, 56);
        if ((uint64_t)(a * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
            return COMPOSITE_FOR_SURE; // divisible by 17
        if ((uint64_t)(a * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
            return COMPOSITE_FOR_SURE; // divisible by 29
        if ((uint64_t)(a * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
            return COMPOSITE_FOR_SURE; // divisible by 43
        if ((uint64_t)(a * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
            return COMPOSITE_FOR_SURE; // divisible by 113
        if ((uint64_t)(a * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
            return COMPOSITE_FOR_SURE; // divisible by 127

        // 2^36-1 is divisible by 3,5,7,19,37,73,109, ...
        a = mpz_mod_mersenne(n, 36);
        if ((uint64_t)(a * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
            return COMPOSITE_FOR_SURE; // divisible by 19
        if ((uint64_t)(a * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
            return COMPOSITE_FOR_SURE; // divisible by 37
        if ((uint64_t)(a * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
            return COMPOSITE_FOR_SURE; // divisible by 73
        if ((uint64_t)(a * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
            return COMPOSITE_FOR_SURE; // divisible by 109

        // 2^44-1 is divisible by 3, 5, 23, 89, .....
        a = mpz_mod_mersenne(n, 44);
        if ((uint64_t)(a * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
            return COMPOSITE_FOR_SURE; // divisible by 23
        if ((uint64_t)(a * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
            return COMPOSITE_FOR_SURE; // divisible by 89

        // 2^23-1 is divisible by 47, .....
        a = mpz_mod_mersenne(n, 23);
        if ((uint64_t)(a * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
            return COMPOSITE_FOR_SURE; // divisible by 47

        // 2^52-1 is divisible by 3, 5, 53, 157, .....
        a = mpz_mod_mersenne(n, 52);
        if ((uint64_t)(a * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
            return COMPOSITE_FOR_SURE; // divisible by 53

        // next small prime to test 59
    }

    // no trivial small factor detected
    return UNDECIDED; // might be prime
}

// integer square root (rounded down)
// assume x < 2^63
static uint64_t uint64_isqrt(uint64_t x)
{
    // Avoid divide by zero
    if (x < 2)
    {
        return x;
    }
    // This code is based on the fact that
    // sqrt(x) == x^1/2 == 2^(log2(x)/2)
    //
    uint64_t log2x = uint64_log_2(x);
    uint64_t log2y = log2x / 2;
    uint64_t y = 1ul << log2y;
    uint64_t y_squared = 1ul << (2 * log2y);
    int64_t sqr_diff = x - y_squared;
    // Perform lerp between powers of four
    y += (sqr_diff / 3) >> log2y;
    // The estimate is probably too low, refine it upward
    y_squared = y * y;
    sqr_diff = x - y_squared;
    y += sqr_diff / (2 * y);
    // The estimate may be too high. If so, refine it downward
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    // The estimate may still be too high
    y -= (-sqr_diff / (2 * y)) + 1;
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    // The oscillating estimate may still be too high
    y -= (-sqr_diff / (2 * y)) + 1;
    y_squared = y * y;
    sqr_diff = x - y_squared;
    if (sqr_diff >= 0)
    {
        return y;
    }
    // The estimate may still be 1 too high
    return sqr_diff < 0 ? y - 1 : y;
}

static bool uint64_is_perfect_square(uint64_t a)
{
    if (0xffedfdfefdecull & (1ull << (a % 48)))
        return false;
    if (0xfdfdfdedfdfcfdecull & (1ull << (a % 64)))
        return false;
    if (0x7bfdb7cfedbafd6cull & (1ull << (a % 63)))
        return false;
    if (0x7dcfeb79ee35ccull & (1ull << (a % 55)))
        return false;
    if (0x8ec196bf5a60dc4ull & (1ull << (a % 61)))
        return false;
    if (0x5d49de7c1846d44ull & (1ull << (a % 59)))
        return false;
    if (0xd228fccfc512cull & (1ull << (a % 53)))
        return false;
    if (0x7bcae4d8ac20ull & (1ull << (a % 47)))
        return false;
    if (0x4a77c5c11acull & (1ull << (a % 43)))
        return false;
    if (0x4c7d4af8c8ull & (1ull << (a % 41)))
        return false;
    if (0x9a1dee164ull & (1ull << (a % 37)))
        return false;
    if (0x6de2b848ull & (1ull << (a % 31)))
        return false;
    if (0xc2edd0cull & (1ull << (a % 29)))
        return false;
    if (0x7acca0ull & (1ull << (a % 23)))
        return false;
    if (0x4f50cull & (1ull << (a % 19)))
        return false;
    if (0x5ce8ull & (1ull << (a % 17)))
        return false;
    if (0x9e4ull & (1ull << (a % 13)))
        return false;

    uint64_t r = uint64_isqrt(a);
    return r * r == a;
}

static bool mpz_is_perfect_square(mpz_t n)
{
    mpz_t ignore;
    mpz_init(ignore);
    uint64_t a = mpz_mod_ui(ignore, n, 64ull * 63ull * 55ull * 61ull * 59ull * 53ull * 47ull * 43ull * 41ull * 37ull);
    uint64_t b = mpz_mod_ui(ignore, n, 31ull * 29ull * 23ull * 19ull * 17ull * 13ull);
    mpz_clear(ignore);

    if (0xffedfdfefdecull & (1ull << (a % 48)))
        return false;
    if (0xfdfdfdedfdfcfdecull & (1ull << (a % 64)))
        return false;
    if (0x7bfdb7cfedbafd6cull & (1ull << (a % 63)))
        return false;
    if (0x7dcfeb79ee35ccull & (1ull << (a % 55)))
        return false;
    if (0x8ec196bf5a60dc4ull & (1ull << (a % 61)))
        return false;
    if (0x5d49de7c1846d44ull & (1ull << (a % 59)))
        return false;
    if (0xd228fccfc512cull & (1ull << (a % 53)))
        return false;
    if (0x7bcae4d8ac20ull & (1ull << (a % 47)))
        return false;
    if (0x4a77c5c11acull & (1ull << (a % 43)))
        return false;
    if (0x4c7d4af8c8ull & (1ull << (a % 41)))
        return false;
    if (0x9a1dee164ull & (1ull << (a % 37)))
        return false;
    if (0x6de2b848ull & (1ull << (b % 31)))
        return false;
    if (0xc2edd0cull & (1ull << (b % 29)))
        return false;
    if (0x7acca0ull & (1ull << (b % 23)))
        return false;
    if (0x4f50cull & (1ull << (b % 19)))
        return false;
    if (0x5ce8ull & (1ull << (b % 17)))
        return false;
    if (0x9e4ull & (1ull << (b % 13)))
        return false;

    // compute the rounded-down integer square root, check if the solution is exact.
    mpz_t u;
    mpz_init(u);
    uint64_t e = mpz_root(u, n, 2); // e = squareroot(n), e non zero if computation is exact.
    mpz_clear(u);
    return e ? true : false;
}

static bool uint64_cipolla(uint64_t n, bool verbose = false)
{
    uint64_t a, d, s = uint64_isqrt(n);
    uint128_t t, u;

    // verify the input number is not a perfect square, where the search of quadratic non-residue would fail
    if (s * s == n)
    {
        return false; // n is a perfect square, and is composite
    }

    // search the smallest non-trivial quadratic residue
    do
    {
        s += 1;
        d = squaremod(s, n);               // d is a quadratic residue
    } while (uint64_is_perfect_square(d)); // d is non-trivial

    // search a quadratic non-residue of the form s^2-d
    while (1)
    {
        s += 1;
        t = (uint128_t)s * s - d;
        a = longlongmod(t, n);
        int j = uint64_jacobi(a, n);
        if (j == 0)
        {
            return false; // composite
        }
        if (j == -1)
        {
            break; // a is a quadratic non-residue
        }
    }

    // Cipolla modular square root finding
    uint64_t e = (n + 1) / 2;
    uint64_t si = 1;
    uint64_t sr = s;
    uint64_t j = 1ul << (uint64_log_2(e) - 1);
    while (j)
    {
        t = (uint128_t)sr * sr + (uint128_t)a * squaremod(si, n);
        si = mulmod(si << 1, sr, n);
        sr = longlongmod(t, n);
        if (e & j)
        {
            t = (uint128_t)s * sr + (uint128_t)a * si;
            u = (uint128_t)s * si + sr;
            sr = longlongmod(t, n);
            si = longlongmod(u, n);
        }
        j = j >> 1;
    }

    // if the outcome (sr, si) is not the square root of d, then n is composite for sure
    if (si)
    {
        return false; // composite for sure
    }
    return (squaremod(sr, n) == d);
}

static bool mpz_cipolla(mod_precompute_t *p, bool verbose = false)
{
    mpz_t a, d, s, t;
    int k;
    unsigned bit;
    bool b = false;
    mpz_inits(a, d, s, t, 0);

    mpz_root(s, p->m, 2);

    // verify the input number is not a perfect square, where the search of quadratic non-residue would fail
    int j = mpz_root(s, p->m, 2);
    if (j != 0)
    {
        // input modulus is a perfect square
        if (verbose)
        {
            printf("Input number is a perfect square\n");
        }
        goto done;
    }

    // search the smallest non-trivial quadratic residue
    do
    {
        mpz_add_ui(s, s, 1);            // s += 1;
        mpz_mul(d, s, s);               // d = s * s
        mpz_mod(d, d, p->m);            // d %= n, d is a quadratic residue
    } while (mpz_is_perfect_square(d)); // d is non-trivial

    // search a quadratic non-residue of the form s^2-d
    while (1)
    {
        mpz_add_ui(s, s, 1); // s += 1
        mpz_mul(t, s, s);    // t = s * s - d
        mpz_sub(t, t, d);
        mpz_mod(a, t, p->m); // a = t % n
        k = mpz_jacobi(a, p->m);
        if (k == 0)
        {
            if (verbose)
            {
                printf("Input number has a factor\n");
            }
            goto done; // composite
        }
        if (k == -1)
        {
            break; // a is a quadratic non-residue
        }
    }

    // gmp_printf("modulus ........ : 0x%Zx\n", p->m);
    // gmp_printf("square d ....... : 0x%Zx\n", d);
    // gmp_printf("non-residue a .. : 0x%Zx\n", a);
    // gmp_printf("s .............. : 0x%Zx\n", s);

    // Cipolla modular square root finding
    //
    mpz_t e, si, sr, u;
    mpz_inits(e, si, sr, u, 0);
    mpz_add_ui(e, p->m, 1);
    mpz_div_2exp(e, e, 1); // e = (n + 1) / 2
    mpz_set_ui(si, 1);
    mpz_set(sr, s);
    bit = mpz_sizeinbase(e, 2) - 1;
    mpz_mod_to_montg(sr, p);
    mpz_mod_to_montg(si, p);
    while (bit--)
    {
        // t = sr * sr + a * si * si;
        mpz_mul(u, sr, sr);
        mpz_mul(t, si, si);
        mpz_mul(t, t, a);
        mpz_add(t, t, u);
        // u = 2 * si * sr
        mpz_mul(u, si, sr);
        mpz_mul_2exp(si, u, 1);
        // si = u % m
        mpz_mod_fast_reduce(si, u, p);
        // sr = t % m
        mpz_set(sr, t);
        mpz_mod_fast_reduce(sr, t, p);

        if (mpz_tstbit(e, bit))
        {
            // t = s * sr + a * si;
            mpz_mul(u, s, sr);
            mpz_mul(t, a, si);
            mpz_add(t, t, u);
            // u = s * si + sr;
            mpz_mul(u, s, si);
            mpz_add(u, u, sr);
            // si = u % m
            // sr = t % m
            if (p->montg)
            {
                // should not use the montgomery reduction here, because this
                // is not after a multiplication of numbers in montgomery form.
                mpz_mod(si, u, p->m);
                mpz_mod(sr, t, p->m);
            }
            else
            {
                mpz_set(si, u);
                mpz_mod_fast_reduce(si, u, p);
                mpz_set(sr, t);
                mpz_mod_fast_reduce(sr, t, p);
            }
        }
    }
    // final reductions, make sure sr and si are < modulus
    mpz_mod_from_montg(sr, t, p);
    mpz_mod_from_montg(si, u, p);
    mpz_mod_slow_reduce(sr, p->m);
    mpz_mod_slow_reduce(si, p->m);

    // if the outcome (sr, si) is not the square root of d, then input number is composite for sure
    // check si == 0
    b = mpz_sgn(si) == 0;
    if (b)
    {
        // compute sr^2
        mpz_mul(t, sr, sr);
        mpz_mod(sr, t, p->m);
        // check sr^2 is the square root of d
        b = mpz_cmp(sr, d) == 0;
        if (verbose)
        {
            if (b)
            {
                printf("Modular square root found\n");
            }
            else
            {
                printf("Modular square root NOT found\n");
            }
        }
    }
    else
    {
        if (verbose)
        {
            printf("Field element is not an integer\n");
        }
    }
    mpz_clears(e, si, sr, u, 0);

done:
    mpz_clears(a, d, s, t, 0);
    return b;
}

static bool uint64_cipolla_primality(uint64_t n, bool verbose = false)
{
    if (n >> 61)
    {
        // the cipolla test might overflow for numbers > 61 bits along
        // the inner additions.
        // Better use another slower deterministic test for numbers <= 60 bits
        // More precise constraint is n < 2^64/6
        mpz_t v;
        mpz_init_set_ui(v, n);
        bool r = mpz_cipolla_primality(v, verbose);
        mpz_clear(v);
        return r;
    }

    if ((n & 1) == 0)
    {
        if (verbose)
        {
            printf("Number is even\n");
        }
        return n == 2; // even
    }

    // sieve small numbers with small factors
    sieve_t sv = uint64_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        if (verbose)
        {
            printf("Number has a small factor\n");
        }
        return false; // composite
    case PRIME_FOR_SURE:
        if (verbose)
        {
            printf("Small number has no small factor\n");
        }
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    bool b = uint64_cipolla(n, verbose);
    if (b != true)
    {
        if (verbose)
        {
            printf("Number is not a Lucas PRP\n");
        }
        return false; // composite
    }
    // prime (BPSW proven to 2^64 and higher, no known countereaxample)
    return true;
}

bool mpz_cipolla_primality(mpz_t n, bool verbose)
{
    if (verbose)
    {
        gmp_printf("Testing a %lu digits number\n", mpz_sizeinbase(n, 10));
    }

    if (mpz_cmp_ui(n, 1ull << 61) < 0)
    {
        // the cipolla test will run in 64 bits calculations
        return uint64_cipolla_primality(mpz_get_ui(n), verbose);
    }

    if (mpz_tstbit(n, 0) == 0)
    {
        if (verbose)
        {
            printf("Number is even\n");
        }
        return false; // even
    }

    // detects small primes, small composites
    // detects smooth composites
    sieve_t sv = mpz_composite_sieve(n);
    switch (sv)
    {
    case COMPOSITE_FOR_SURE:
        if (verbose)
        {
            printf("Number has a small factor\n");
        }
        return false; // composite
    case PRIME_FOR_SURE:
        if (verbose)
        {
            printf("Small number has no small factor\n");
        }
        return true; // prime
    case UNDECIDED:
    default:
        break;
    }

    mod_precompute_t *pcpt = mpz_mod_precompute(n, verbose);
    bool r = mpz_cipolla(pcpt, verbose);
    mpz_mod_uncompute(pcpt);

    if (verbose && r == false)
    {
        printf("Number is composite\n");
    }
    if (verbose && r == true)
    {
        printf("Number passed all tests and is unlikely a composite one\n");
    }

    return r;
}

// ------------------------------------------------------------------------------
// Simple foolguard unit tests
// ------------------------------------------------------------------------------

void cipolla_primality_self_test(void)
{
    uint64_t a, b, m;
    uint128_t aa;

    // ---------------------------------------------------------------------------------
    printf("Modular arithmetic sanity check ...\n");

    m = 0x8765432187654321ull;
    a = 0xffffffffffffffffull;
    b = longmod(0, a, m);
    assert(b == 8690466094656961758ull);
    aa = a;
    b = longlongmod(aa, m);
    assert(b == 8690466094656961758ull);

    b = longmod(1, a, m);
    assert(b == 7624654210261333660ull);
    aa = ((uint128_t)1 << 64) + a;
    b = longlongmod(aa, m);
    assert(b == 7624654210261333660ull);

    // ---------------------------------------------------------------------------------
    printf("Fast reduction mod 2^b-1 ...\n");
    mpz_t ma, mb;
    mpz_inits(ma, mb, 0);
    mpz_set_ui(ma, 11);
    for (unsigned b = 1; b < 10; b++)
    {
        mpz_mul(ma, ma, ma);
    }
    for (unsigned b = 1; b < 64; b++)
    {
        uint64_t mm = (1ull << b) - 1;
        uint64_t r = mpz_mod_mersenne(ma, b);
        assert(r % mm == mpz_mod_ui(mb, ma, mm));
    }

    // ---------------------------------------------------------------------------------
    printf("Sieve ...\n");
    assert(uint64_composite_sieve(101) == PRIME_FOR_SURE);
    assert(uint64_composite_sieve(1661) == COMPOSITE_FOR_SURE);
    assert(uint64_composite_sieve(281474976710677ull) == UNDECIDED);

    // 2^127 - 1 (a prime)
    mpz_init_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 127);
    mpz_sub_ui(ma, ma, 1);
    assert(mpz_composite_sieve(ma) == UNDECIDED);
    // 2^127 + 1
    mpz_add_ui(ma, ma, 2);
    assert(mpz_composite_sieve(ma) == COMPOSITE_FOR_SURE);

    // ---------------------------------------------------------------------------------
    printf("Slow modular reduction\n");
    mpz_t x;
    mpz_init(x);
    mpz_set_ui(ma, 0x1);
    mpz_mul_2exp(ma, ma, 127);
    mpz_add_ui(ma, ma, 0x1d);
    // verify 2*(modulus +1) == 2
    mpz_add_ui(mb, ma, 0x1);
    mpz_add(x, mb, mb);
    mpz_mod_slow_reduce(x, ma);
    assert(mpz_cmp_ui(x, 2) == 0);

    // ---------------------------------------------------------------------------------
    printf("Fast modular reduction\n");
    // verify generic modulus
    mpz_set_ui(ma, 0x1cc);
    mod_precompute_t *p = mpz_mod_precompute(ma);
    assert(p->special_case == false);
    mpz_set_ui(mb, 0x2000ull * 0x1cc);
    mpz_mod_fast_reduce(mb, x, p);
    mpz_mod_slow_reduce(mb, ma);
    assert(mpz_get_ui(mb) == 0);
    // clears temp structure
    mpz_mod_uncompute(p);

    // verify proth modulus
    mpz_t mx, mtt;
    mpz_init(mx);
    mpz_init(mtt);
    mpz_set_ui(ma, 0xcdefabcdcdefabcdull);
    mpz_mul_2exp(ma, ma, 64);
    mpz_add_ui(ma, ma, 1);
    p = mpz_mod_precompute(ma);
    assert(p->special_case == true);
    assert(p->montg == true);
    assert(p->proth == true);
    assert(p->n2 == 64);
    assert(p->n == 128);
    mpz_set_ui(ma, 0xabcdef01abcdef01ull);
    mpz_mul(ma, ma, ma);
    mpz_set_ui(mb, 0x1234567812345678ull);
    mpz_mul(mb, mb, mb);
    mpz_mul(x, ma, mb);
    mpz_mod(x, x, p->m);
    mpz_mod_to_montg(ma, p);
    mpz_mod_to_montg(mb, p);
    mpz_mul(mx, ma, mb);
    mpz_mod_fast_reduce(mx, mtt, p);
    mpz_mod_from_montg(mx, mtt, p);
    mpz_mod_slow_reduce(mx, p->m);
    assert(mpz_cmp(x, mx) == 0);
    // clears temp structure
    mpz_clear(mx);
    mpz_mod_uncompute(p);

    // verify (modulus + 0x17)^2 == 17*17
    mpz_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 127);
    p = mpz_mod_precompute(ma);
    assert(p->special_case == false);
    mpz_add_ui(mb, ma, 17);
    mpz_mul(x, mb, mb);
    mpz_mod_fast_reduce(x, mtt, p);
    mpz_mod_slow_reduce(x, ma);
    assert(mpz_get_ui(x) == 17 * 17);

    // verify (modulus + 0x19)^4 == 17*17*17*17
    mpz_add_ui(mb, ma, 19);
    mpz_mul(x, mb, mb);
    mpz_mul(x, x, x);
    mpz_mod_fast_reduce(x, mtt, p);
    mpz_mod_slow_reduce(x, ma);
    assert(mpz_get_ui(x) == 19 * 19 * 19 * 19);

    // clears temp structure
    mpz_mod_uncompute(p);
    p = 0;
    mpz_clear(x);

    // ---------------------------------------------------------------------------------
    printf("Small primes (uint64)\n");
    assert(uint64_cipolla_primality(16777259ull) == true);
    assert(uint64_cipolla_primality(281474976710677ull) == true);

    // ---------------------------------------------------------------------------------
    printf("Small composites (uint64)\n");
    assert(uint64_cipolla_primality(16777265ull) == false);
    assert(uint64_cipolla_primality(281474976710683ull) == false);

    // ---------------------------------------------------------------------------------
    printf("Small primes (uint64 sanity check)\n");

    bool res = uint64_cipolla_primality(31, false);
    assert(res == true);

    res = uint64_cipolla_primality(65537, false);
    assert(res == true);

    // ---------------------------------------------------------------------------------
    printf("Small primes (mpz sanity check)\n");
    mpz_t ms, me, mtmp;
    mpz_inits(ms, me, mtmp, 0);

    // quick tests mod 31
    mpz_set_ui(ma, 31);
    p = mpz_mod_precompute(ma);

    res = mpz_cipolla_primality(ma);
    assert(res == true);

    mpz_mod_uncompute(p);

    // quick tests mod 2^89+29
    mpz_set_ui(ma, 1);
    mpz_mul_2exp(ma, ma, 89);
    mpz_add_ui(ma, ma, 29);
    p = mpz_mod_precompute(ma);

    res = mpz_cipolla_primality(ma);
    assert(res == true);

    mpz_mod_uncompute(p);

    // quick tests mod 21*2^128+1
    mpz_set_ui(ma, 21);
    mpz_mul_2exp(ma, ma, 128);
    mpz_add_ui(ma, ma, 1);
    p = mpz_mod_precompute(ma);
    assert(p->proth == true);

    res = mpz_cipolla_primality(ma);
    assert(res == true);

    mpz_mod_uncompute(p);
    mpz_clears(ms, me, mtmp, 0);

    // ---------------------------------------------------------------------------------
    uint32_t smallq[] = {
        1,   1,   1,   3,   1,   5,   3,   3,   1,  9,   7,   5,   3,   17,  27,  3,   1,   29,  3,   21,  7,   17,
        15,  9,   43,  35,  15,  29,  3,   11,  3,  11,  15,  17,  25,  53,  31,  9,   7,   23,  15,  27,  15,  29,
        7,   59,  15,  5,   21,  69,  55,  21,  21, 5,   159, 3,   81,  9,   69,  131, 33,  15,  135, 29,  13,  131,
        9,   3,   33,  29,  25,  11,  15,  29,  37, 33,  15,  11,  7,   23,  13,  17,  9,   75,  3,   171, 27,  39,
        7,   29,  133, 59,  25,  105, 129, 9,   61, 105, 7,   255, 277, 81,  267, 81,  111, 39,  99,  39,  33,  147,
        27,  51,  25,  281, 43,  71,  33,  29,  25, 9,   451, 41,  277, 165, 67,  27,  7,   29,  51,  17,  169, 39,
        67,  27,  27,  33,  85,  155, 87,  155, 37, 5,   217, 5,   175, 27,  85,  51,  91,  69,  147, 45,  253, 95,
        27,  15,  45,  69,  97,  299, 7,   107, 19, 21,  117, 141, 85,  83,  87,  147, 49,  129, 105, 77,  7,   9,
        427, 75,  87,  309, 15,  165, 49,  215, 27, 159, 205, 303, 57,  35,  129, 5,   133, 65,  27,  35,  21,  107,
        15,  101, 235, 351, 67,  15,  7,   581, 33, 203, 375, 47,  33,  71,  57,  75,  7,   251, 423, 129, 163, 185,
        217, 81,  49,  189, 735, 119, 735, 483, 3,  249, 67,  105, 357, 431, 43,  81,  25,  249, 67,  29,  115, 261,
        69,  59,  133, 315, 337, 63,  81,  119, 25, 65,  421, 39,  79,  95,  297, 155, 73,  435, 223, 0};

    printf("Medium primes (mpz)\n");
    for (int j = 0; smallq[j]; j++)
    {
        // primes (2^j + q) must be catched
        mpz_set_ui(ma, 1);
        mpz_mul_2exp(ma, ma, j);
        mpz_add_ui(ma, ma, smallq[j]);
        bool res = mpz_cipolla_primality(ma, false);
        if (!res)
        {
            gmp_printf("prime 0x%Zx\n", ma);
            printf("error at %d %d\n", (int)j, (int)smallq[j]);
            assert(res == true);
        }
    }

    // ---------------------------------------------------------------------------------
    uint32_t smallp[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37,  41,  43,
                         47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 0};
    printf("Medium composites (mpz)\n");
    for (int j = 0; smallp[j]; j++)
    {
        // Composites p * (2^127-1) must be catched
        mpz_set_ui(ma, 1);
        mpz_mul_2exp(ma, ma, 127);
        mpz_sub_ui(ma, ma, 1);
        mpz_mul_ui(ma, ma, smallp[j]);
        bool res = mpz_cipolla_primality(ma);
        if (res)
        {
            printf("error at %d %d\n", (int)j, (int)smallp[j]);
            assert(res == false);
        }
    }

    // ---------------------------------------------------------------------------------
    printf("Large primes (mpz)\n");

    // large proth prime 333*2^448+1
    mpz_set_ui(ma, 333);
    mpz_mul_2exp(ma, ma, 448);
    mpz_add_ui(ma, ma, 1);
    assert(mpz_cipolla_primality(ma) == true);

    // 11111...6442446...11111 (1001-digits) The smallest zeroless titanic palindromic prime
    // https://t5k.org/curios/page.php?number_id=3797
    char titanic[1002];
    memset(titanic, '1', 1001);
    titanic[497] = '6';
    titanic[498] = '4';
    titanic[499] = '4';
    titanic[500] = '2';
    titanic[501] = '4';
    titanic[502] = '4';
    titanic[503] = '6';
    titanic[1001] = 0;
    mpz_set_str(ma, titanic, 10);
    // gmp_printf("%Zd\n", ma);
    assert(mpz_cipolla_primality(ma) == true);

    // 11111...0...3781 (1001-digits) The smallest Cyclops titanic prime.
    // https://t5k.org/curios/page.php?number_id=18967
    memset(titanic, '1', 1001);
    titanic[500] = '0';
    titanic[997] = '3';
    titanic[998] = '7';
    titanic[999] = '8';
    titanic[1000] = '1';
    titanic[1001] = 0;
    mpz_set_str(mb, titanic, 10);
    assert(mpz_cipolla_primality(mb) == true);

    // ---------------------------------------------------------------------------------
    printf("Large composites (mpz)\n");
    // a semiprime out of the 2 previous tests, no small factors.
    mpz_mul(ma, mb, ma);
    assert(mpz_cipolla_primality(ma) == false);

    // a large square
    mpz_mul(ma, mb, mb);
    assert(mpz_cipolla_primality(ma) == false);

    // a large cube
    mpz_mul(ma, ma, mb);
    assert(mpz_cipolla_primality(ma) == false);

    mpz_clears(ma, mb, 0);
}
