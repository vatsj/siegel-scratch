"""Interval-arithmetic certification of the explicit Green argument.

Re-checks, with outward-rounded interval arithmetic (mpmath.iv), everything
that green_explicit.py / tail.py / certify.py did in floating point:

  1. each interval [d_lo, d_hi] of results/certificate.json: the Lemma 1
     criterion (exact integers) and margin > 0, i.e.
        2 sum_{p odd inert} K_p log p  >  (Lemma 2 bound)
     for every class-number-one d in the interval;
  2. the tail d >= 10^40: the linear bound margin >= 0.18 log d - 16 > 0;
  3. the gap (163, d0]: an exact integer split-prime check.

Floats are used only to choose evaluation points (the Abel grid).  Every
inequality is then evaluated in intervals and holds whatever those points are.

Usage: python3 rigorous.py
"""

import json
import math
import os

import mpmath
import numpy as np
from mpmath import iv

from certify import CERT, small_d_check
from green_explicit import P_ENUM, Params, K_exact, _primes_upto, injective
from tail import D_TAIL, T_DROP, tail_params

iv.prec = 160
mpmath.mp.prec = 160

PRIMES = [int(p) for p in _primes_upto(P_ENUM)]
LOG_P = [iv.log(p) for p in PRIMES]
THETA_ENUM = sum(LOG_P[1:], LOG_P[0])  # interval containing theta(P_ENUM)
GRID = 40000


def lo(x):
    return mpmath.mpf(x.a)


def hi(x):
    return mpmath.mpf(x.b)


def X(v):
    """Exact interval for an int, a float, or an mpf."""
    return iv.mpf(v)


# ---------------------------------------------------------------- Lemma 2

def arch_upper(P, d_hi):
    """Interval whose upper end bounds the Lemma 2 bound / (R^2 S) at d_hi."""
    R, S, J = X(P.R), X(P.S), X(P.J)
    a, b = P.N1 // 2, P.N2 // 2
    rho2 = X(a * a) + X(d_hi) * X(b * b)
    lr = iv.log(rho2) / 2
    val = lr * (R + S - 2) / R + iv.log(X(P.N1 * P.N2 * P.J)) / R
    log_R = iv.log(R)
    if lo(log_R) > hi(lr):  # R > rho for certain
        val += 2 * (J - 1) / R * (log_R - lr - 1 + iv.exp(lr - log_R))
    elif hi(log_R) > lo(lr):
        raise ValueError("cannot decide R > rho")
    return val


# ---------------------------------------------------------------- Lemma 3

def _G_point(Y, R, S):
    """Interval for G(Y) at an exact point Y (mpf)."""
    S_ = X(S)
    A = (S_ - X(0.5)) / 2
    cap = S - 1
    if Y <= 0:
        return X(0)
    Yi = X(Y)
    if Y <= cap:
        return A * Yi ** 2 - Yi ** 3 / 6
    capi = X(cap)
    G_cap = A * capi ** 2 - capi ** 3 / 6
    return G_cap + (Yi - capi) * S_ * (S_ - 1) / 2


def f_interval(t, R, S):
    """Interval containing f(t) = t G((R-1)/t - 1).  G is nondecreasing, so we
    evaluate it at the two ends of the interval for Y."""
    ti = X(t)
    Y = (X(R) - 1) / ti - 1
    G_lo = lo(_G_point(lo(Y), R, S))
    G_hi = hi(_G_point(hi(Y), R, S))
    return ti * iv.mpf([G_lo, G_hi])


def theta_lower(t):
    """Rosser-Schoenfeld: theta(t) > t (1 - 1/(2 log t)) for t >= 563."""
    ti = X(t)
    assert t >= 563
    return lo(ti * (1 - 1 / (2 * iv.log(ti))))


def weighted_prime_sum_lower(R, S):
    """Rigorous lower bound (mpf) for sum_{3 <= p < R} K_p log p."""
    Xmax = R - 1
    total = X(0)
    for p, lp in zip(PRIMES[1:], LOG_P[1:]):
        if p > Xmax:
            break
        total += X(K_exact(p, R, S)) * lp
    total = lo(total)
    if Xmax <= P_ENUM:
        return total
    # Abel summation on (P_ENUM, R-1] (f(R-1) = 0, f decreasing, theta
    # nondecreasing):
    #   sum >= -f(P_ENUM) theta(P_ENUM) + sum_k theta(t_k) (f(t_k) - f(t_k+1))
    pts = [mpmath.mpf(P_ENUM)]
    pts += [mpmath.mpf(float(x)) for x in np.geomspace(P_ENUM, float(Xmax), GRID)[1:-1]]
    pts.append(mpmath.mpf(Xmax))
    pts = sorted(set(pts))
    fs = [f_interval(t, R, S) for t in pts]
    acc = -hi(fs[0] * THETA_ENUM)
    for k in range(len(pts) - 1):
        diff = lo(fs[k] - fs[k + 1])
        if diff > 0:
            acc += lo(X(theta_lower(pts[k])) * X(diff))
    return total + acc


def split_upper(log_x, log_d):
    """Interval version of reduction.split_upper (log_x, log_d intervals)."""
    log_b = iv.log(X(2)) + (log_x - log_d) / 2
    if hi(log_b) < 0:
        return X(0)
    log_b = iv.mpf([max(lo(log_b), mpmath.mpf(0)), hi(log_b)])
    inv_sqrt_d = iv.exp(-log_d / 2)
    return 4 * log_x * (iv.pi ** 2 / 6 * inv_sqrt_d ** 2 + iv.pi * inv_sqrt_d * (1 + log_b))


def nonarch_lower(P, d_lo):
    """Rigorous lower bound (mpf) for 2 sum_{p odd inert} K_p log p / (R^2 S),
    for every class-number-one d >= d_lo."""
    R, S = P.R, P.S
    total = X(weighted_prime_sum_lower(R, S))
    log_d = iv.log(X(d_lo))
    k_max = X(S) * X(R) * (X(R) - 1) / 2  # K_p <= k_max / p
    total -= k_max * split_upper(iv.log(X(R - 1)), log_d)
    if d_lo < R:
        total -= k_max * log_d / X(d_lo)
    return lo(2 * total / (X(R) ** 2 * X(S)))


def check_interval(rec):
    P = Params(*(int(rec[k]) for k in ("N1", "N2", "R", "S", "J")))
    assert injective(P), "Lemma 1 criterion fails"
    assert rec["d_lo"] >= 164 and rec["d_lo"] <= rec["d_hi"]
    m = nonarch_lower(P, rec["d_lo"]) - hi(arch_upper(P, rec["d_hi"]))
    return m


# ------------------------------------------------------------------ tail

def check_tail():
    """Interval check of the linear tail bound derived in NOTES.md.

    For d >= 10^40 (L = log d) the derivation gives
        margin >= slope * L + intercept - (1e-38 L + 1e-8)
    with the slope/intercept below.  The last bracket absorbs every o(1) term:
    log(1 - 4/d), log(1 - 64/R), (R-1)^2/R^2, 3L/R, and the 1e-9 corrections in
    the bounds for log rho (each is < 1e-38 L or < 1e-9 for d >= 10^40).
    """
    T = X(T_DROP)
    eta = 1 / (2 * iv.log(T))
    theta_T = sum((LOG_P[i] for i in range(len(PRIMES)) if PRIMES[i] <= T_DROP), X(0))
    boundary = (theta_T - (1 - eta) * T) / (2 * T)  # >= f(T)(theta(T)-(1-eta)T)/(R^2 S)
    head = (1 + iv.log(T / 64)) / 2  # >= int_0^T f / (R^2 S), using S <= R/64
    # main >= (1/2) log S - 1/2 >= (1/2) L - (1/2) log 256 - 1/2
    nonarch_slope = (1 - eta)  # 2 (1 - eta) * (1/2)
    nonarch_const = 2 * (1 - eta) * (-iv.log(X(256)) / 2 - X(0.5) - head) - 2 * boundary
    c_rho = iv.log(X(0.138125)) / 2  # rho^2 <= 0.138125 d^{3/2} (1 + 1e-9)
    arch_slope = X(65) / 64 * X(0.75) + X(1) / 64  # (65/64) lr + (1/16)(L/4)
    arch_const = X(65) / 64 * c_rho
    slope = nonarch_slope - arch_slope
    intercept = nonarch_const - arch_const
    L0 = iv.log(X(D_TAIL))
    s_lo, c_lo = lo(slope), lo(intercept)
    # need: slope L + intercept - 1e-38 L - 1e-8 >= 0.18 L - 16 for all L >= L0
    ok = (s_lo - mpmath.mpf("0.18") - mpmath.mpf("1e-38") > 0
          and (s_lo - mpmath.mpf("0.18") - mpmath.mpf("1e-38")) * lo(L0)
          + c_lo + 16 - mpmath.mpf("1e-8") >= 0
          and mpmath.mpf("0.18") * lo(L0) - 16 > 0)
    # injectivity for all d >= 10^40 (NOTES.md): N1 (N2 - 33.01) > 1.0157 R,
    # which follows from N1 N2 >= 1.1 R and N2 >= 10^10 / 2.
    ok_inj = lo(X(1.1) * (1 - X(33.01) / X(5 * 10 ** 9))) > mpmath.mpf("1.0157")
    return ok and ok_inj, s_lo, c_lo


# ------------------------------------------------------------------ main

def main(path=CERT):
    with open(path) as fh:
        cert = json.load(fh)
    ivs = sorted(cert["intervals"], key=lambda r: r["d_lo"])
    for a, b in zip(ivs, ivs[1:]):
        assert a["d_hi"] >= b["d_lo"], ("gap", a, b)
    assert ivs[-1]["d_hi"] >= D_TAIL, "certificate must reach the tail"
    worst = None
    for rec in ivs:
        m = check_interval(rec)
        status = "ok" if m > 0 else "FAIL"
        print(f"[{rec['d_lo']:.6g}, {rec['d_hi']:.6g}]  rigorous margin >= "
              f"{mpmath.nstr(m, 6)}  (float {rec['margin']:+.6f})  {status}", flush=True)
        assert m > 0, rec
        worst = m if worst is None else min(worst, m)
    ok, s, c = check_tail()
    print(f"tail d >= 1e40: slope >= {mpmath.nstr(s, 8)}, intercept >= "
          f"{mpmath.nstr(c, 8)} -> {'ok' if ok else 'FAIL'}")
    assert ok
    d0 = ivs[0]["d_lo"]
    p = small_d_check(d0)
    print(f"small d: every prime d = 3 mod 8 in (163, {d0:.0f}] has a split prime "
          f"p <= {p} below d/4")
    print(f"CERTIFIED: no imaginary quadratic field with d > 163 has class number one "
          f"(min interval margin {mpmath.nstr(worst, 4)}).")


if __name__ == "__main__":
    main()
