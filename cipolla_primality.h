#pragma once

// -----------------------------------------------------------------------
// Cubic primality test
//
// mpz_cipolla_primality():
//    true: composite for sure
//    false: might be prime
//
// cipolla_primality_self_test()
//    simplified unit tests to detect a possible compiler/platform issue.
//    assert when fail (this should not happen).
// -----------------------------------------------------------------------

#include "gmp.h"
#include <stdbool.h>

bool mpz_cipolla_primality(mpz_t v, bool verbose = false);
void cipolla_primality_self_test(void);
