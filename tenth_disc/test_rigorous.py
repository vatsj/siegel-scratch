"""Checks for rigorous.py (run: python3 -m unittest -v test_rigorous)."""

import json
import unittest

import numpy as np

from certify import CERT
from green_explicit import P_ENUM, _primes_upto


class TestRigorous(unittest.TestCase):
    def test_rosser_schoenfeld_on_enumerated_range(self):
        # theta(t) > t (1 - 1/(2 log t)) for t >= 563.  theta is a right-
        # continuous step function, so on [563, P_ENUM] the worst points are
        # t = 563 itself and the left limits just before each prime > 563.
        p = _primes_upto(P_ENUM)
        theta = np.cumsum(np.log(p.astype(float)))
        bound = lambda t: t * (1 - 1 / (2 * np.log(t)))
        i563 = int(np.searchsorted(p, 563))
        self.assertEqual(p[i563], 563)
        self.assertGreater(theta[i563], bound(563.0))
        t = p[i563 + 1:].astype(float)
        self.assertTrue(np.all(theta[i563:-1] > bound(t)))

    def test_tail_certificate(self):
        from rigorous import check_tail

        ok, slope, intercept = check_tail()
        self.assertTrue(ok)
        self.assertGreater(slope, 0.18)
        self.assertGreater(intercept, -16)

    def test_one_certificate_interval(self):
        from rigorous import check_interval

        with open(CERT) as fh:
            recs = json.load(fh)["intervals"]
        top = max(recs, key=lambda r: r["d_hi"])
        self.assertGreater(check_interval(top), 0)


if __name__ == "__main__":
    unittest.main()
