# CipollaPrimalityTest

WORK IN PROGRESS !!!

This is test code to verify a cipolla primality test, based on linear recurrences 

So far, no counterexample (pseudoprime) has been found.

 

# The maths

Cipolla algorithm, when p is prime, is able to retrieve a modular square root

https://en.wikipedia.org/wiki/Cipolla%27s_algorithm

The converse is : If a modular square root cannot be retrieved, then p is composite for sure

https://math.stackexchange.com/questions/3551965/can-cipolla-algorithm-be-used-as-a-primality-test

There might be an infinity of pseudoprimes for this test, but none has been found.

No false positive, no false negative found < 2^48

# See also

See the following self-similar primality tests https://github.com/Boutoukoat

- SimplePrimalityTest
- CipollaPrimalityTest
- QuadraticPrimalityTest
- CubicPrimalityTest

# Cipolla utility based of GMP library for large integers

```
$ make

$ make check

$ ./cipolla 2^127-1 2^127+0 2^127+1
2^127-1 might be prime, time=       0.151 msecs.
2^127+0 is composite for sure, time=       0.000 msecs.
2^127+1 is composite for sure, time=       0.000 msecs.

$ ./cipolla 0x988a04da39838a3757afef4ae6ed84b092aa0ee673067e52140862e5d27af3adfd1d65489e91b068df21f5de5e78fe4a8deb967201c7944b0a0eabc31bb0b824d3cb6293156c0c84bc48072952f08711da7a8786050335f82ec0bba57adf9c22aad36ba2f4919a3ccd8a4717799d90ffc82189f5425a3026de65b4c7e11e9beb
0x988a04da39838a3757afef4ae6ed84b092aa0ee673067e52140862e5d27af3adfd1d65489e91b068df21f5de5e78fe4a8deb967201c7944b0a0eabc31bb0b824d3cb6293156c0c84bc48072952f08711da7a8786050335f82ec0bba57adf9c22aad36ba2f4919a3ccd8a4717799d90ffc82189f5425a3026de65b4c7e11e9beb might be prime, time=       6.194 msecs.

$ ./cipolla 2^11213-1
2^11213-1 might be prime, time=    3147.283 msecs.

```

# Complete user's guide :

Later.

# Limits

GMP and memory limits apply.








