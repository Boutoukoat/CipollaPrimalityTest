#pragma once

// -----------------------------------------------------------------------
// Parse an expression (string) and convert to a mpz_t
//
// The expression can be very long and complex (large numbers, brackets ...).
// Limits come from GMP.
//
// e.g. "2*(12*12+1)-1" --> 289
// -----------------------------------------------------------------------

#include "gmp.h"
void mpz_expression_parse(mpz_t n, char *str);
