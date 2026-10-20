# siegel-scratch

Scratch work for Terence Tao's
[experimental AI-math challenge](https://terrytao.wordpress.com/2026/10/20/an-experimental-ai-math-challenge-using-tenth-discriminant-bounds-to-explore-arguments-ruling-out-siegel-zeroes/)
(20 Oct 2026): for each new argument ruling out Landau–Siegel zeroes, extract an
explicit upper bound on a hypothetical *tenth* imaginary quadratic field
Q(√−d) of class number one, without Heegner's methods, and optimise it.
Background: Tao's [sketch of Green's argument](https://terrytao.wordpress.com/2026/10/19/a-sketch-of-the-proof-of-non-existence-of-siegel-zeroes/).

## Step 1: the reduction "explicit no-Siegel bound ⇒ bound on d"

Every such argument ends in (or can be pushed to) a bound of the form

    (H)   Σ_{p ≤ x, χ(p) = −1} log p / p  ≤  α log x + β log q + γ log log x + C

(Pescador–Barretto: α = ½, β = ¼, γ = 1, C = O(1)). The remaining step is to turn
(H) into a bound on d. [`tenth_disc/reduction.py`](tenth_disc/reduction.py) does
this explicitly.

**Observation.** If h(−d) = 1 with d > 163, then d is prime, d ≡ 3 (mod 8), q = d,
and every *split* prime is a norm p = (a² + d b²)/4 with b ≠ 0. So split primes are
≥ d/4 **and** very sparse beyond that:

    Σ_{p ≤ x split} log p / p  ≤  4 log x · ( (π²/6)/d + (π/√d)(1 + log(2√(x/d))) ).

Together with Rosser–Schoenfeld's explicit Mertens bound
Σ_{p≤x} log p/p > log x − 1.3326 − 1/(2 log x), the inert sum is ≈ log x − 1.33
up to x ≈ exp(d^{1/4}), not just up to x = d/4. Plugging a huge x into (H) gives
a bound on d that is *polynomial* in the constant C, instead of exponential.

| (α, β, γ) = (½, ¼, 1) | Silberman's reduction (x = q/16 only) | split-prime reduction |
|---|---|---|
| C = 0    | — | d ≤ 10^5.4 |
| C = 10   | q ≤ 10^44.6 (with C_M = 2, as in his comment) | d ≤ 10^6.4 |
| C = 100  | — | d ≤ 10^8.1 |
| C = 1000 | — | d ≤ 10^10.0 |

Bounds of 10^6 to 10^10 are well within direct computation of class numbers, so
**the whole difficulty now sits in making C explicit** (and checking the range of
x on which (H) is proved; `NoSiegelBound.min_power` handles bounds valid only for
x ≥ d^A).

Caveats: the threshold is found by a numerical scan over log d and log x, not
yet a certified interval computation. Also, this reduction is elementary and
classical in spirit (Heilbronn–Linfoot-style use of the principal form), so
there is no claim of novelty.

## Next steps

- Make Green's bound (2) fully explicit: the O(ε) terms in Lemma 2 (Hadamard),
  Lemma 3 (the k_{i,j} basis count), the O(1) losses in Lemma 1 (resultant degree
  count), and the optimisation over ε, N, R, S, J. That gives (H) with an explicit C.
- Do the same for the Banks and Kerr write-ups and the QRH papers, and compare.
- Make the threshold rigorous (interval arithmetic or a monotonicity proof).

## Running

    cd tenth_disc
    python3 reduction.py                     # table above
    python3 -m unittest -v test_reduction    # checks bounds vs exact sums for d = 11..163

Standard library only (Python ≥ 3.9).

*AI disclosure: the code and write-up were produced with Claude Code (Anthropic).*
