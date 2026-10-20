"""Sanity checks for reduction.py (run: python3 -m unittest discover tenth_disc)."""

import math
import unittest

from reduction import (
    NoSiegelBound,
    best_margin,
    inert_lower,
    silberman_threshold_log_q,
    threshold_log_d,
)

HEEGNER_ODD = (11, 19, 43, 67, 163)  # class number one, d = 3 mod 8, d prime


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if sieve[i]]


def chi(d, p):
    """Kronecker symbol (-d / p) for d = 3 mod 8 prime."""
    if p == 2:
        return -1
    if p == d:
        return 0
    return 1 if pow(-d % p, (p - 1) // 2, p) == 1 else -1


class TestInertLowerBound(unittest.TestCase):
    """For genuine h = 1 fields the bound must hold for every x."""

    X = 10 ** 6

    def test_against_exact_sums(self):
        ps = primes_upto(self.X)
        for d in HEEGNER_ODD:
            total, idx = 0.0, 0
            for x in sorted({3, 10, 40, 100, 1000, 10 ** 4, 10 ** 5, self.X}):
                while idx < len(ps) and ps[idx] <= x:
                    p = ps[idx]
                    if chi(d, p) == -1:
                        total += math.log(p) / p
                    idx += 1
                lower = inert_lower(math.log(x), math.log(d))
                self.assertLessEqual(lower, total, f"d={d}, x={x}")

    def test_primes_below_d_over_4_are_inert(self):
        for d in HEEGNER_ODD:
            for p in primes_upto(d // 4):
                self.assertEqual(chi(d, p), -1, f"d={d}, p={p}")


class TestReductions(unittest.TestCase):
    def test_silberman_reproduces_10_to_45(self):
        L = silberman_threshold_log_q(0.5, 0.25, C_LS=10.0, C_M=2.0)
        self.assertAlmostEqual(L / math.log(10), 44.6, delta=0.1)

    def test_split_prime_reduction_beats_silberman(self):
        h = NoSiegelBound(alpha=0.5, beta=0.25, gamma=1.0, C=10.0)
        L = threshold_log_d(h)
        self.assertIsNotNone(L)
        self.assertLess(L / math.log(10), 7.0)
        self.assertGreater(best_margin(h, L + 1.0)[0], 0.0)

    def test_trivial_bound_excludes_nothing(self):
        # alpha = 1 is just Mertens: it can never contradict h = 1.
        h = NoSiegelBound(alpha=1.0, beta=0.0, gamma=0.0, C=0.0)
        self.assertIsNone(threshold_log_d(h))


if __name__ == "__main__":
    unittest.main()
