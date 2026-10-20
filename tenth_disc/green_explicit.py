"""Green's polynomial-method argument, made explicit for the tenth discriminant.

Notation follows Tao's sketch
(terrytao.wordpress.com/2026/10/19/a-sketch-of-the-proof-of-non-existence-of-siegel-zeroes/),
with D = -d.  Parameters: a grid size N1 x N2, degrees (R, S), derivative
order J.

  * Lambda = integer polynomials P(z1, z2) with deg_{z1} < R, deg_{z2} < S.
  * E = {(a, conj a) : a = m + n sqrt(D)}, with m in a centred window of
    N1 consecutive integers and n in a centred window of N2.
  * T P = (d^j/dz1^j P (e))_{j < J, e in E}  in  Z[sqrt D]^{J x |E|}.
  * M = Gram determinant of T(Lambda) (real inner product on C = R^2).
    M is a positive integer when T is injective.

Assume h(-d) = 1 with d > 163, so d is prime, d = 3 mod 8, every prime
p < d/4 is inert, and split primes are norms (a^2 + d b^2)/4 (see
reduction.py).  Then

    sum_{p odd inert} 2 K_p log p  <=  log M  <=  ArchUpper,          (*)

where K_p = sum_{i<R, j<S} min(floor(i/p), S-1-j) (Lemma 3, exact count).
If (*) fails for some admissible parameters, then h(-d) = 1 is impossible.

Lemma 1 (injectivity), general form.  T is injective whenever

    J (N1 N2 - S + 1)  >  2 (R-1)(S-1) + max(N1, N2) (R + S - 2).

Proof sketch (no irreducibility assumption).  Let P != 0 vanish to
horizontal order >= J on E, and write P = A(z2) prod_i F_i^{e_i}, with
F_i distinct irreducible factors of z1-degree r_i >= 1 and z2-degree s_i.
  - A has degree <= S-1 and the z2-coordinates of E are distinct, so A
    vanishes at <= S-1 points of E.  Let E' be the other points.
  - At e = (z0, c) in E', sum_i e_i mu_i(e) >= J, where mu_i(e) is the
    order of z0 as a root of F_i(., c).
  - Res_{z1}(F_i, dF_i/dz1) is a nonzero polynomial in z2 of degree
    <= (2 r_i - 1) s_i.  It is divisible by (z2 - c)^{mu_i(e) - 1} for
    each e (multilinearity of the Sylvester determinant, as in Tao's
    sketch).  Hence sum_e (mu_i(e) - 1)_+ <= (2 r_i - 1) s_i.
  - #{e : mu_i(e) >= 1} <= (r_i + s_i) max(N1, N2): apply Schwartz-Zippel
    to the nonzero polynomial F_i(m + n sqrt D, m - n sqrt D) in (m, n).
  Summing e_i times these bounds over i, with sum e_i r_i <= R-1 and
  sum e_i s_i <= S-1:
      J |E'| <= 2 (R-1)(S-1) + max(N1, N2) (R + S - 2).

Lemma 2 (Archimedean).  Hadamard with the monomial basis, |alpha| <= rho
on E, and falling factorial a^(j) <= a^j give

    log M <= R S log(|E| J) + R S (R + S - 2) log rho
             + 2 S (J-1) (R log(R/rho) - R + rho)          [if R > rho]

with rho^2 = (N1//2)^2 + d (N2//2)^2.

Lemma 3 (non-Archimedean).  For odd inert p, z1^p - z2 lies in
p Z[sqrt D] on E, and every z1-derivative of order >= 1 of z1^p - z2 is
divisible by p.  The basis z1^{i - p k} (z1^p - z2)^k z2^j, with
k = min(floor(i/p), S-1-j), is unitriangular against the monomial basis.
Hence v_p(M) >= 2 K_p.

Prime sums use exact enumeration for p <= P_ENUM.  Beyond that they use Abel
summation against a decreasing minorant f(t) <= K_t, together with
Rosser-Schoenfeld (1962): t (1 - 1/(2 log t)) < theta(t) for t >= 563.

Floating point is used throughout, so the results are "numerically
explicit", not certified by interval arithmetic.
"""

import math
from dataclasses import dataclass

import numpy as np

from reduction import split_upper

P_ENUM = 2_000_000


def _primes_upto(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0]


_PRIMES = _primes_upto(P_ENUM)
_PRIMES_ODD = _PRIMES[1:]
_LOG_PRIMES_ODD = np.log(_PRIMES_ODD.astype(float))
_THETA_ENUM = float(np.log(_PRIMES.astype(float)).sum())  # theta(P_ENUM)


@dataclass(frozen=True)
class Params:
    N1: int
    N2: int
    R: int
    S: int
    J: int


def injective(P):
    """Explicit Lemma 1 criterion."""
    lhs = P.J * (P.N1 * P.N2 - P.S + 1)
    rhs = 2 * (P.R - 1) * (P.S - 1) + max(P.N1, P.N2) * (P.R + P.S - 2)
    return lhs > rhs


def max_R(N1, N2, S, J):
    """Largest R satisfying the Lemma 1 criterion (or 0)."""
    Nm = max(N1, N2)
    num = J * (N1 * N2 - S + 1) + 2 * (S - 1) - Nm * (S - 2)
    den = 2 * (S - 1) + Nm
    if num <= 0:
        return 0
    R = -(-num // den) - 1  # largest R with den * R < num
    return max(R, 0)


def log_rho(P, d):
    """log rho, rho^2 = (N1//2)^2 + d (N2//2)^2, computed without overflow."""
    a, b = P.N1 // 2, P.N2 // 2
    terms = []
    if a:
        terms.append(2 * math.log(a))
    if b:
        terms.append(math.log(d) + 2 * math.log(b))
    hi = max(terms)
    return 0.5 * (hi + math.log(sum(math.exp(t - hi) for t in terms)))


def arch_norm(P, d):
    """Lemma 2 upper bound for log M, divided by R^2 S."""
    R, S, J = P.R, P.S, P.J
    lr = log_rho(P, d)
    log_R = math.log(R)
    val = lr * ((R + S - 2) / R) + math.log(P.N1 * P.N2 * J) / R
    if log_R > lr:
        val += 2 * ((J - 1) / R) * (log_R - lr - 1 + math.exp(lr - log_R))
    return val


def arch_upper(P, d):
    """Lemma 2 upper bound for log M (unnormalised; small parameters only)."""
    return arch_norm(P, d) * (P.R ** 2 * P.S)


def K_exact(p, R, S):
    """sum_{i<R, j<S} min(floor(i/p), S-1-j), exactly (Python ints)."""

    def g(q):
        return q * (S - 1) - q * (q - 1) // 2 if q <= S - 1 else S * (S - 1) // 2

    def gsum(n):  # sum_{q<n} g(q)
        m = min(n, S)
        s = (S - 1) * m * (m - 1) // 2 - m * (m - 1) * (m - 2) // 6
        if n > S:
            s += (n - S) * S * (S - 1) // 2
        return s

    Q, r0 = divmod(R, p)
    return p * gsum(Q) + r0 * g(Q)


def _K_vec(p, R, S):
    """K_exact vectorised over a numpy array of primes, in float64 (R < 2^53)."""
    Q = (R // p).astype(float)
    r0 = (R % p).astype(float)
    S = float(S)

    def g(q):
        return np.where(q <= S - 1, q * (S - 1) - q * (q - 1) / 2, S * (S - 1) / 2)

    m = np.minimum(Q, S)
    gs = (S - 1) * m * (m - 1) / 2 - m * (m - 1) * (m - 2) / 6
    gs += np.where(Q > S, (Q - S) * S * (S - 1) / 2, 0.0)
    return p * gs + r0 * g(Q)


def _f_norm(s, R, S):
    """f(sR) / (R^2 S), where f(t) = t G((R-1)/t - 1) <= K_t is decreasing.

    G(Y) = int_0^Y g~, with g~(y) = y (S - 1/2) - y^2/2 on [0, S-1] (the
    continuous extension of q -> sum_j min(q, S-1-j)), constant afterwards.
    Written so that no intermediate quantity overflows for huge R.
    """
    R, S = float(R), float(S)
    Y = (1.0 - 1.0 / R) / s - 1.0
    A = (S - 0.5) / 2.0
    cap = S - 1.0
    sY, yR = s * Y, Y / R
    region1 = sY * yR * (A / S) - sY * yR * (Y / S) / 6.0
    region2 = (s * cap) * (cap / R) * ((A - cap / 6.0) / S) + (sY - s * cap) * ((S - 1.0) / (2.0 * R))
    return np.where(Y <= 0, 0.0, np.where(Y <= cap, region1, region2))


def _f_minorant(t, R, S):
    """Unnormalised f(t) (small parameters only; used by the tests)."""
    return _f_norm(np.asarray(t, dtype=float) / R, R, S) * (float(R) ** 2 * S)


def weighted_prime_sum_norm(R, S, grid=40000):
    """Lower bound for sum_{3 <= p < R} K_p log p, divided by R^2 S."""
    X = R - 1
    small = _PRIMES_ODD[_PRIMES_ODD <= min(X, P_ENUM)]
    if R < 2 ** 53:
        K_small = _K_vec(small, R, S) / (float(R) ** 2 * S)
    else:  # beyond exact float/int64 range: use the minorant f(p) <= K_p
        K_small = _f_norm(small.astype(float) / float(R), R, S)
    total = float((K_small * _LOG_PRIMES_ODD[: len(small)]).sum())
    if X <= P_ENUM:
        return total
    # Abel summation on (P_ENUM, X], with f(X) = 0 and theta nondecreasing:
    #   sum f(p) log p >= -f(P_ENUM) theta(P_ENUM)
    #                     + sum_k theta_lo(t_k) (f(t_k) - f(t_{k+1})).
    Rf = float(R)
    s = np.geomspace(P_ENUM / Rf, X / Rf, grid)
    f = _f_norm(s, R, S)
    log_t = math.log(Rf) + np.log(s)
    theta_lo_over_R = s * (1.0 - 1.0 / (2.0 * log_t))
    total += -float(f[0]) * _THETA_ENUM + Rf * float(
        (theta_lo_over_R[:-1] * (f[:-1] - f[1:])).sum()
    )
    return total


def nonarch_norm(P, d_lo):
    """Lower bound for sum_{p odd inert} 2 K_p log p, divided by R^2 S, valid
    for every class-number-one d >= d_lo (d_lo > 163)."""
    R, S = P.R, P.S
    total = weighted_prime_sum_norm(R, S)
    log_d = math.log(d_lo)
    k_max = (R - 1) / (2 * R)  # K_p <= k_max R^2 S / p
    total -= k_max * split_upper(math.log(R - 1), log_d)  # split primes (>= d/4)
    if d_lo < R:
        total -= k_max * log_d / d_lo  # the ramified prime p = d
    return 2.0 * total


def margin(P, d_lo, d_hi=None):
    """(nonarch lower bound - arch upper bound) / (R^2 S).  Positive means
    every class-number-one d in [d_lo, d_hi] is impossible."""
    if d_hi is None:
        d_hi = d_lo
    if not injective(P):
        return -math.inf
    return nonarch_norm(P, d_lo) - arch_norm(P, d_hi)


def _scale(n, x):
    """floor(n * x) for a big int n and a moderate float x >= 0, exactly enough."""
    return (n * int(x * 2 ** 40)) >> 40


def params_from(log_N1, log_aspect, log_eps, kappa, d, R_cap=None):
    """Build Params: N2 ~ N1 * aspect / sqrt(d), S ~ eps * R, J ~ kappa * S,
    R maximal for Lemma 1 (solved by a short fixed-point iteration)."""
    N1 = max(2, int(math.exp(log_N1)))
    N2 = max(1, int(math.exp(log_N1 + log_aspect - 0.5 * math.log(d))))
    eps = math.exp(log_eps)
    R = N1 * N2
    if R_cap is not None:
        R = min(R, R_cap)
    for _ in range(30):
        S = max(2, _scale(R, eps))
        J = max(2, _scale(S, kappa))
        R_new = max_R(N1, N2, S, J)
        if R_cap is not None:
            R_new = min(R_new, R_cap)
        if R_new <= 2 or R_new == R:
            R = R_new
            break
        R = R_new
    if R <= 2:
        return None
    S = max(2, _scale(R, eps))
    J = max(2, _scale(S, kappa))
    R = max_R(N1, N2, S, J)
    if R_cap is not None:
        R = min(R, R_cap)
    if R <= S:
        return None
    return Params(N1, N2, R, S, J)


def optimise(d_lo, d_hi=None, start=None, iters=60, R_cap=None):
    """Coordinate search over (log N1, log aspect, log eps, kappa)."""
    if d_hi is None:
        d_hi = d_lo
    x = list(start) if start else [0.5 * math.log(d_hi) + 4.0, 0.0, -3.0, 1.0]
    steps = [1.0, 0.5, 0.5, 0.3]

    def score(v):
        P = params_from(*v, d_hi, R_cap=R_cap)
        return (margin(P, d_lo, d_hi), P) if P else (-math.inf, None)

    best, bestP = score(x)
    for _ in range(iters):
        improved = False
        for k in range(4):
            for sgn in (1, -1):
                y = list(x)
                y[k] += sgn * steps[k]
                if k == 3 and y[k] <= 0:
                    continue
                s, P = score(y)
                if s > best:
                    best, bestP, x, improved = s, P, y, True
        if not improved:
            steps = [s / 2 for s in steps]
            if max(steps) < 1e-3:
                break
    return best, bestP, x
