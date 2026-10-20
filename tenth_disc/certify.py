"""Cover [d0, D_MAX] by intervals on which Green's explicit inequality (*)
fails, and close the gap (163, d0] by an elementary check.

Usage:
    python3 certify.py            # search, then write results/certificate.json
    python3 certify.py --verify   # re-check an existing certificate

Output intervals [d_lo, d_hi] each come with explicit parameters
(N1, N2, R, S, J) for which margin(P, d_lo, d_hi) > 0, i.e. no d in the
interval can have class number one (see green_explicit.py).
"""

import json
import math
import os
import sys

from green_explicit import Params, injective, margin, optimise

D_MAX = 10.0 ** 200
R_CAP = 10 ** 300  # keeps every float conversion finite
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "results", "certificate.json")


def _best(d_lo, d_hi, warm):
    """Optimised margin for [d_lo, d_hi]: warm start plus a few cold starts."""
    starts = [warm] if warm else []
    L = math.log(d_hi)
    starts += [[0.5 * L + k, 0.0, e, 1.5] for k in (8, 16, 30) for e in (-3.0, -5.0)]
    best = (-math.inf, None, None)
    for x0 in starts:
        r = optimise(d_lo, d_hi, start=x0, iters=40, R_cap=R_CAP)
        if r[0] > best[0]:
            best = r
        if warm and best[0] > 0.05:
            break  # warm start already good enough
    return best


def search(d_floor=10.0 ** 5):
    """Walk down from D_MAX; each step takes the widest interval it can."""
    intervals = []
    d_hi, warm = D_MAX, None
    while d_hi > d_floor:
        lo, hi = math.log(max(d_floor, 2.0)), math.log(d_hi)
        # largest log-width whose optimised margin is still positive
        width, found = min(hi - lo, 40.0), None
        while width > 1e-3:
            d_lo = math.exp(hi - width)
            m, P, x = _best(d_lo, d_hi, warm)
            if m > 0:
                found = (d_lo, m, P, x)
                break
            width /= 2
        if not found:
            break
        d_lo, m, P, x = found
        intervals.append({"d_lo": d_lo, "d_hi": d_hi, "margin": m,
                          "N1": str(P.N1), "N2": str(P.N2), "R": str(P.R),
                          "S": str(P.S), "J": str(P.J)})
        print(f"[{d_lo:.4g}, {d_hi:.4g}]  margin {m:+.4f}", flush=True)
        d_hi, warm = d_lo, x
    return intervals


def small_d_check(d0):
    """Every prime d = 3 mod 8 with 163 < d <= d0 has a prime p < d/4 with
    (-d/p) = +1, hence h(-d) > 1.  Returns the largest p ever needed."""
    import numpy as np

    from green_explicit import _primes_upto

    ps = [int(p) for p in _primes_upto(2000)[1:]]  # odd primes
    worst = 0
    for d in (int(x) for x in _primes_upto(int(d0)) if x > 163 and x % 8 == 3):
        for p in ps:
            if 4 * p >= d:
                raise AssertionError(f"no split prime below d/4 found for d={d}")
            if pow(-d % p, (p - 1) // 2, p) == 1:
                worst = max(worst, p)
                break
        else:
            raise AssertionError(f"search bound too small for d={d}")
    return worst


def verify(path=CERT):
    with open(path) as fh:
        cert = json.load(fh)
    ivs = sorted(cert["intervals"], key=lambda r: r["d_lo"])
    for a, b in zip(ivs, ivs[1:]):
        assert a["d_hi"] >= b["d_lo"], ("gap", a, b)
    assert ivs[-1]["d_hi"] >= cert["D_MAX"]
    worst = math.inf
    for r in ivs:
        P = Params(*(int(r[k]) for k in ("N1", "N2", "R", "S", "J")))
        assert injective(P), r
        m = margin(P, r["d_lo"], r["d_hi"])
        assert m > 0, (r, m)
        worst = min(worst, m)
    d0 = ivs[0]["d_lo"]
    p = small_d_check(d0)
    print(f"OK: {len(ivs)} intervals cover [{d0:.6g}, {cert['D_MAX']:.3g}], "
          f"min margin {worst:.4f}; all d in (163, {d0:.6g}] killed by a split "
          f"prime p <= {p}.")


if __name__ == "__main__":
    if "--verify" in sys.argv:
        verify()
        sys.exit()
    ivs = search()
    os.makedirs(os.path.dirname(CERT), exist_ok=True)
    with open(CERT, "w") as fh:
        json.dump({"D_MAX": D_MAX, "intervals": ivs}, fh, indent=1)
    verify()
