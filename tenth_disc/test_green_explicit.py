"""Brute-force checks of the explicit Lemmas 1-3 on tiny parameters.

Run: python3 -m unittest -v test_green_explicit
"""

import math
import unittest
from itertools import product

import numpy as np

from green_explicit import (
    Params,
    _f_minorant,
    _K_vec,
    K_exact,
    arch_upper,
    injective,
    margin,
    max_R,
    optimise,
)


def legendre_neg_d(d, p):
    return 1 if pow(-d % p, (p - 1) // 2, p) == 1 else -1


def grid(N):
    lo = -(N // 2)
    return range(lo, lo + N)


def mul(x, y, D):
    """(a + b sqrt D)(c + e sqrt D) as integer pairs."""
    return (x[0] * y[0] + D * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def power(x, k, D):
    r = (1, 0)
    for _ in range(k):
        r = mul(r, x, D)
    return r


def T_rows(P, d):
    """Rows T(z1^a z2^b) as integer vectors (re, im-coefficient pairs)."""
    D = -d
    pts = [((m, n), (m, -n)) for m in grid(P.N1) for n in grid(P.N2)]
    rows = []
    for a in range(P.R):
        for b in range(P.S):
            row = []
            for j in range(P.J):
                ff = math.perm(a, j) if j <= a else 0
                for z1, z2 in pts:
                    if ff == 0:
                        row += [0, 0]
                        continue
                    v = mul(power(z1, a - j, D), power(z2, b, D), D)
                    row += [ff * v[0], ff * v[1]]
            rows.append(row)
    return rows


def gram_det(rows, d):
    """det of Gram matrix with <x, y> = sum x_re y_re + d x_im y_im."""
    n = len(rows)
    w = [1, d] * (len(rows[0]) // 2)
    G = [[sum(wi * x * y for wi, x, y in zip(w, r1, r2)) for r2 in rows] for r1 in rows]
    # Bareiss fraction-free elimination
    sign, prev = 1, 1
    for k in range(n - 1):
        if G[k][k] == 0:
            for i in range(k + 1, n):
                if G[i][k] != 0:
                    G[k], G[i] = G[i], G[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                G[i][j] = (G[i][j] * G[k][k] - G[i][k] * G[k][j]) // prev
        prev = G[k][k]
    return sign * G[n - 1][n - 1]


def vp(M, p):
    v = 0
    while M % p == 0:
        M //= p
        v += 1
    return v


CASES = [
    # (d, Params) -- small enough for exact determinants.  The first four
    # satisfy the Lemma 1 criterion; the last does not.
    (7, Params(N1=5, N2=4, R=13, S=2, J=5)),
    (11, Params(N1=6, N2=3, R=9, S=3, J=6)),
    (19, Params(N1=4, N2=4, R=10, S=3, J=6)),
    (43, Params(N1=6, N2=4, R=14, S=2, J=5)),
    (67, Params(N1=3, N2=3, R=8, S=3, J=3)),
]


class TestLemmasBruteForce(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dets = {(d, P): gram_det(T_rows(P, d), d) for d, P in CASES}

    def test_injective_cases_have_nonzero_M(self):
        found = False
        for (d, P), M in self.dets.items():
            if injective(P):
                found = True
                self.assertGreater(M, 0, (d, P))
        self.assertTrue(found, "no case satisfies the Lemma 1 criterion")

    def test_archimedean_bound(self):
        for (d, P), M in self.dets.items():
            if M > 0:
                self.assertLessEqual(math.log(M), arch_upper(P, d) + 1e-9, (d, P))

    def test_nonarchimedean_bound(self):
        checked = 0
        for (d, P), M in self.dets.items():
            if M <= 0:
                continue
            for p in (3, 5, 7, 11, 13):
                if p < P.R and d % p and legendre_neg_d(d, p) == -1:
                    K = K_exact(p, P.R, P.S)
                    self.assertGreaterEqual(vp(M, p), 2 * K, (d, P, p))
                    checked += K > 0
        self.assertGreater(checked, 0)


class TestCountingFormulas(unittest.TestCase):
    def test_K_exact_matches_double_sum(self):
        for R, S, p in product((1, 7, 20, 33), (1, 2, 5, 9), (2, 3, 5, 7, 31)):
            brute = sum(min(i // p, S - 1 - j) for i in range(R) for j in range(S))
            self.assertEqual(K_exact(p, R, S), brute, (R, S, p))
            self.assertAlmostEqual(
                float(_K_vec(np.array([p]), R, S)[0]), brute, delta=1e-6 * max(1, brute)
            )

    def test_minorant_below_K_and_decreasing(self):
        for R, S in ((50, 7), (400, 30), (10 ** 5, 900)):
            ps = np.arange(2, R + 5)
            f = _f_minorant(ps.astype(float), R, S)
            K = _K_vec(ps, R, S)
            self.assertTrue(np.all(f <= K * (1 + 1e-12) + 1e-9), (R, S))
            self.assertTrue(np.all(np.diff(f) <= 1e-9 * np.abs(f[:-1]) + 1e-9), (R, S))

    def test_max_R_is_maximal(self):
        for N1, N2, S, J in ((40, 7, 5, 6), (300, 20, 40, 40), (1000, 3, 9, 30)):
            R = max_R(N1, N2, S, J)
            if R > 0:
                self.assertTrue(injective(Params(N1, N2, R, S, J)))
            self.assertFalse(injective(Params(N1, N2, R + 1, S, J)))


class TestUnconditionalInequality(unittest.TestCase):
    """(*) with the *true* inert primes is a theorem for every d, so the
    inert-weighted sum may never exceed the Archimedean bound."""

    def test_real_splitting_data(self):
        import random

        from green_explicit import P_ENUM, _LOG_PRIMES_ODD, _PRIMES_ODD

        rng = random.Random(1)
        worst, n = 0.0, 0
        for d in (163, 1019, 10007, 99991, 1000003, 3, 7, 43):
            for _ in range(25):
                N1, N2 = rng.randint(10, 3000), rng.randint(1, 300)
                S, J = rng.randint(2, 2000), rng.randint(2, 4000)
                R = max_R(N1, N2, S, J)
                if R <= S or R > P_ENUM:
                    continue
                P = Params(N1, N2, R, S, J)
                ps = _PRIMES_ODD[_PRIMES_ODD < R]
                inert = np.array([d % int(p) != 0 and legendre_neg_d(d, int(p)) == -1
                                  for p in ps])
                K = _K_vec(ps, R, S)
                lhs = 2 * float((K * _LOG_PRIMES_ODD[: len(ps)] * inert).sum())
                worst = max(worst, lhs / arch_upper(P, d))
                n += 1
        self.assertGreater(n, 50)
        self.assertLess(worst, 1.0)


class TestSanity(unittest.TestCase):
    def test_cannot_exclude_genuine_heegner_numbers(self):
        # h(-d) = 1 really holds here, so no admissible parameters may "win".
        # (The large-d inputs -- d prime = 3 mod 8, p < d/4 inert -- are true
        # for these d too.)
        for d in (11, 19, 43, 67, 163):
            best, _, _ = optimise(d)
            self.assertLessEqual(best, 0.0, d)


if __name__ == "__main__":
    unittest.main()
