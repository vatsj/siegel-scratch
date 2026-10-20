"""Checks for tail.py (run: python3 -m unittest -v test_tail)."""

import math
import unittest

import numpy as np

from green_explicit import _PRIMES, injective, margin
from tail import (
    D_TAIL,
    THETA_T,
    T_DROP,
    I_closed,
    linear_lower,
    tail_margin,
    tail_params,
)

SAMPLES = [int(10 ** e) for e in (40, 41.5, 47, 60, 88, 123, 150, 199)] + [
    10 ** 40 + 3, 7 * 10 ** 52 + 19, 10 ** 250 + 11
]


class TestTail(unittest.TestCase):
    def test_theta_constant_is_an_upper_bound(self):
        p = _PRIMES[_PRIMES <= T_DROP].astype(float)
        theta = math.fsum(np.log(p))
        self.assertLessEqual(theta, THETA_T)
        self.assertLess(THETA_T - theta, 0.01)

    def test_I_closed_matches_definition(self):
        # (1/S) int_1^oo G(u-1) u^-3 du by direct quadrature (moderate S).
        for S in (7, 50, 900):
            A = (S - 0.5) / 2
            cap = S - 1
            G_cap = A * cap ** 2 - cap ** 3 / 6
            B = S * (S - 1) / 2
            u = np.geomspace(1.0, float(S), 400001)
            G = A * (u - 1) ** 2 - (u - 1) ** 3 / 6
            region1 = np.trapezoid(G / u ** 3, u)
            region2 = G_cap / (2 * S ** 2) + B / (2 * S)
            self.assertAlmostEqual(I_closed(S), (region1 + region2) / S, places=6)

    def test_I_closed_symbolic(self):
        try:
            import sympy as sp
        except ImportError:
            self.skipTest("sympy not installed")
        S, u = sp.symbols("S u", positive=True)
        A = (S - sp.Rational(1, 2)) / 2
        cap = S - 1
        G1 = A * (u - 1) ** 2 - (u - 1) ** 3 / 6
        Gcap = A * cap ** 2 - cap ** 3 / 6
        B = S * (S - 1) / 2
        I = (sp.integrate(G1 / u ** 3, (u, 1, S))
             + sp.integrate((Gcap + (u - S) * B) / u ** 3, (u, S, sp.oo))) / S
        target = sp.log(S) / 2 - sp.Rational(1, 2) + (sp.log(S) + 2) / (4 * S)
        self.assertEqual(sp.simplify(I - target), 0)

    def test_tail_positive_and_above_linear_bound(self):
        for d in SAMPLES:
            self.assertGreaterEqual(d, D_TAIL)
            P = tail_params(d)
            self.assertTrue(injective(P), d)
            self.assertLess(4 * (P.R - 1), d)  # all odd p < R are < d/4, hence inert
            m = tail_margin(d)
            self.assertGreaterEqual(m, linear_lower(math.log(d)), d)
            self.assertGreater(m, 0.0, d)

    def test_full_engine_agrees(self):
        # The numerical engine uses sharper estimates, so it must do at least
        # as well as the closed form with the same parameters.
        for d in SAMPLES:
            if d > 10 ** 200:
                continue
            P = tail_params(d)
            self.assertGreaterEqual(margin(P, d) + 1e-9, tail_margin(d), d)


if __name__ == "__main__":
    unittest.main()
