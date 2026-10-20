"""Closed-form version of the explicit Green inequality for very large d.

For d >= D_TAIL we use the regime R <= d/4, where every odd prime p < R is
inert (if h(-d) = 1), so neither split primes nor the ramified prime occur.
Everything is then a closed-form expression, and the normalised margin is at
least a linear function of log d with positive slope.  This handles the
infinite tail that a finite certificate cannot reach.

Parameters (all explicit integers):
    R = floor(d/4),  S = floor(R/64),  J = 2 S,
    N2 = ceil(d^{1/4} / 2),  N1 = ceil(1.1 R / N2).

Non-Archimedean side.  Drop the primes p <= T := 10^6 (they contribute >= 0).
Write theta(t) = sum_{p<=t} log p and f(t) = t G((R-1)/t - 1) <= K_t, which
is decreasing with f(R-1) = 0.  Abel summation gives

    sum_{T<p<R} f(p) log p = -f(T) theta(T) + int_T^{R-1} theta(t) (-f'(t)) dt.

For t >= T, Rosser-Schoenfeld gives theta(t) >= (1 - eta) t with
eta = 1/(2 log T).  Hence

    sum >= (1 - eta) (int_0^{R-1} f - int_0^T f) - f(T) (theta(T) - (1 - eta) T).

Here int_0^{R-1} f = (R-1)^2 S I(S), with I(S) in closed form below.  Also
f(t) <= min(R S^2/2, R^2 S/(2t)), which gives
int_0^T f <= (R^2 S/2) (1 + log(T S/R)) and f(T) <= R^2 S/(2T).
"""

import math

from green_explicit import Params, arch_norm, injective

T_DROP = 10 ** 6
THETA_T = 998484.18  # upper bound for theta(10^6) = 998484.17503 (see test_tail.py)
D_TAIL = 10 ** 40


def I_closed(S):
    """(1/S) int_1^oo G(u - 1) u^{-3} du = (log S)/2 - 1/2 + (log S + 2)/(4 S).

    Exact identity (verified symbolically with sympy; see test_tail.py).
    """
    lS = math.log(S)
    return 0.5 * lS - 0.5 + (lS + 2.0) / (4.0 * S)


def tail_params(d):
    d = int(d)
    R = d // 4
    S = R // 64
    J = 2 * S
    N2 = -(-math.isqrt(math.isqrt(d)) // 2) + 1  # >= d^{1/4}/2
    N1 = -(-11 * R // (10 * N2))  # ceil(1.1 R / N2)
    return Params(N1, N2, R, S, J)


def tail_margin(d):
    """Closed-form lower bound for the normalised margin (R <= d/4 regime)."""
    P = tail_params(d)
    R, S = P.R, P.S
    eta = 1.0 / (2.0 * math.log(T_DROP))
    main = ((R - 1) / R) ** 2 * I_closed(S)
    head = 0.5 * (1.0 + math.log(T_DROP * S / R))
    boundary = (THETA_T - (1.0 - eta) * T_DROP) / (2.0 * T_DROP)
    nonarch = 2.0 * ((1.0 - eta) * (main - head) - boundary)
    return nonarch - arch_norm(P, d)


def linear_lower(log_d):
    """Simple linear bound: tail_margin(d) >= 0.18 log d - 16 for d >= D_TAIL
    (see NOTES.md for the derivation)."""
    return 0.18 * log_d - 16.0


if __name__ == "__main__":
    for e in (12, 20, 30, 37, 40, 60, 100, 200, 300):
        d = 10 ** e
        P = tail_params(d)
        print(f"d=1e{e}: injective={injective(P)} closed-form margin={tail_margin(d):+.4f}"
              f"  linear={linear_lower(math.log(d)):+.4f}")
