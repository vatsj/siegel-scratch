"""From an explicit "no Landau-Siegel zero" bound to a tenth-discriminant bound.

Setting.  Suppose Q(sqrt(-d)) has class number one with d > 163.  By genus
theory d is then a prime with d = 3 mod 8, the conductor of chi = (-d/.) is
q = d, and the only ramified prime is d itself.

Every split prime p is the norm of a principal ideal, i.e.

    p = (a^2 + d b^2) / 4   with a, b integers, b != 0,

so split primes are >= d/4 and very sparse (unconditionally, given h = 1).
Hence the "inert mass"

    I(x) := sum_{p <= x, chi(p) = -1} log p / p

is essentially the full Mertens sum log x - E up to x ~ exp(d^{1/4}), and
not merely up to x = d/4 (which is all that a naive argument uses).

Any explicit Siegel-zero-free upper bound of the shape

    I(x) <= alpha log x + beta log q + gamma log log x + C     (H)

with alpha < 1 then contradicts h = 1 once d is large enough.  This module
computes that threshold.  Everything here is elementary and fully explicit;
the hard part of the challenge -- proving (H) with explicit constants -- is
*not* addressed here.

Explicit inputs used:
  * Rosser-Schoenfeld (1962), (3.21): for x > 1,
        sum_{p <= x} log p / p > log x - E - 1/(2 log x),
    with E = 1.33258227573322...
  * sum_{a in Z} 1/(a^2 + c^2) <= 1/c^2 + pi/c   (compare with the integral).

All quantities are parametrised by log x and log d so that astronomically
large x, d can be handled with plain floats.
"""

import math
from dataclasses import dataclass

# E = -gamma - sum_p log p / (p (p - 1)), Rosser-Schoenfeld notation (sign flipped).
E_MERTENS = 1.332582275733220881765828776071027748


def mertens_lower(log_x):
    """Rosser-Schoenfeld lower bound for sum_{p <= x} log p / p (x > 1)."""
    return log_x - E_MERTENS - 1.0 / (2.0 * log_x)


def split_upper(log_x, log_d):
    """Upper bound for sum_{p <= x split} log p / p, assuming h(-d) = 1.

    Each split p <= x equals (a^2 + d b^2)/4 with 1 <= b <= 2 sqrt(x/d), and
    log p / p <= 4 log x / (a^2 + d b^2).  Summing over all a in Z and those b:

        <= 4 log x * sum_b ( 1/(d b^2) + pi/(sqrt(d) b) )
        <= 4 log x * ( (pi^2/6)/d + pi/sqrt(d) * (1 + log B) ),

    with B = 2 sqrt(x/d).
    """
    log_b = math.log(2.0) + 0.5 * (log_x - log_d)
    if log_b < 0.0:
        return 0.0
    inv_sqrt_d = math.exp(-0.5 * log_d)
    return 4.0 * log_x * (
        (math.pi ** 2 / 6.0) * inv_sqrt_d ** 2 + math.pi * inv_sqrt_d * (1.0 + log_b)
    )


def inert_lower(log_x, log_d):
    """Lower bound for I(x) = sum_{p <= x, chi(p) = -1} log p / p if h(-d) = 1."""
    ramified = log_d * math.exp(-log_d) if log_x >= log_d else 0.0
    return mertens_lower(log_x) - ramified - split_upper(log_x, log_d)


@dataclass(frozen=True)
class NoSiegelBound:
    """A hypothesised explicit bound (H), valid for x >= d^min_power."""

    alpha: float
    beta: float
    gamma: float
    C: float
    min_power: float = 0.0

    def rhs(self, log_x, log_d):
        return (
            self.alpha * log_x
            + self.beta * log_d
            + self.gamma * math.log(log_x)
            + self.C
        )


def best_margin(bound, log_d, samples=400):
    """max over admissible x of inert_lower(x) - rhs(x).

    A positive value means d is excluded: h(-d) = 1 would violate (H) at the
    maximising x.  We scan log x geometrically from max(d^min_power, e^1.01)
    up to roughly exp(d^{1/2}), beyond which the split-prime bound is useless.
    """
    lo = max(bound.min_power * log_d, 1.01)
    hi = max(lo * 1.0001, math.exp(0.5 * log_d) + log_d)
    best_val, best_log_x = -math.inf, lo
    for k in range(samples + 1):
        log_x = lo * (hi / lo) ** (k / samples)
        val = inert_lower(log_x, log_d) - bound.rhs(log_x, log_d)
        if val > best_val:
            best_val, best_log_x = val, log_x
    return best_val, best_log_x


def threshold_log_d(bound, log_d_min=math.log(164), log_d_max=600.0, step=0.005):
    """Smallest L such that every d with log d in [L, log_d_max] is excluded.

    Scans log d downward from log_d_max, first on a coarse grid and then on a
    grid of the given step inside the coarse cell where exclusion first fails.
    Returns None if d near log_d_max is still not excluded (i.e. (H) is too
    weak).  Numerical: exclusion is only checked at grid points.
    """

    def excluded(log_d):
        return best_margin(bound, log_d)[0] > 0.0

    if not excluded(log_d_max):
        return None
    coarse = 0.5
    hi = log_d_max
    while hi - coarse >= log_d_min and excluded(hi - coarse):
        hi -= coarse
    lo = max(hi - coarse, log_d_min)
    # hi is excluded; refine downward from hi towards lo.
    while hi - step >= lo and excluded(hi - step):
        hi -= step
    return hi


def silberman_threshold_log_q(alpha, beta, C_LS, C_M, extra=1.0):
    """Lior Silberman's reduction (blog comment, 20 Oct 2026), x = q/16 only:

        (1 - alpha - beta) log q <= 2 log log q + C_LS + 2 C_M
                                    + (1 - alpha) log 16 + extra

    (extra = 1 is the "+1" he later noted was lost).  Returns the largest
    log q satisfying it, found by fixed-point iteration.
    """
    k = 1.0 - alpha - beta
    if k <= 0:
        return None
    const = C_LS + 2 * C_M + (1 - alpha) * math.log(16) + extra
    L = 10.0
    for _ in range(200):
        L = (2 * math.log(L) + const) / k
    return L


if __name__ == "__main__":
    # Silberman's sample parameters, compared under both reductions.
    a, b, C = 0.5, 0.25, 10.0
    L_lior = silberman_threshold_log_q(a, b, C_LS=C, C_M=2.0)
    print(f"Silberman reduction, alpha={a}, beta={b}, C_LS={C}, C_M=2:")
    print(f"  q <= e^{L_lior:.1f} = 10^{L_lior / math.log(10):.1f}")
    for C in (0.0, 10.0, 100.0, 1000.0):
        h = NoSiegelBound(alpha=a, beta=b, gamma=1.0, C=C)
        L = threshold_log_d(h)
        _, log_x = best_margin(h, L)
        print(
            f"Split-prime reduction, C={C:>6}: d <= 10^{L / math.log(10):.2f}"
            f"  (optimal x ~ e^{log_x:.3g})"
        )
