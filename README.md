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

## Step 2: Green's argument made explicit → d < 515,960

Details and proofs are in [`tenth_disc/NOTES.md`](tenth_disc/NOTES.md).
Rather than passing through (H), I apply Green's lattice argument directly to
a hypothetical tenth discriminant. Lemmas 1–3 of Tao's sketch are all made
explicit; Lemma 1 is proved in general, without the "P irreducible" cheat.
The result is a finite inequality

    2 Σ_{p odd inert} K_p log p ≤ log M ≤ (Hadamard bound)

that h(−d) = 1 forces. Under h(−d) = 1, all p < d/4 are inert and split primes
are sparse, so the left side is large and the inequality fails once d is large.

| range | excluded by |
|---|---|
| (163, 515,960] | elementary: some prime p ≤ 43 below d/4 splits |
| [515,960, 10^200] | 18 intervals with explicit parameters ([`results/certificate.json`](tenth_disc/results/certificate.json)) |
| [10^40, ∞) | closed form: margin ≥ 0.18 log d − 16 ([`tail.py`](tenth_disc/tail.py)) |

So Green's argument recovers the class number one theorem with threshold
≈ 5·10^5, compared with Oesterlé's e^55 ≈ 10^24. Everything is re-checked in
**outward-rounded interval arithmetic** (`rigorous.py`): all 18 margins are
certified positive (smallest 0.00172), and so is the tail bound.

**Checks:**
- brute-force exact Gram determinants confirm Lemmas 1–3 on small cases (the
  v_p bound is often sharp);
- with real splitting data, the unconditional inequality never fails across
  200+ random parameter sets (ratio ≤ 0.52, matching Green's "½");
- genuine Heegner numbers (d ≤ 163) are never excluded;
- the tail closed form agrees with the full engine.

**Caveats:**
- Lemma 1 (general case) has been checked against Green's note (his Lemma
  3.1 and A.1): it is the same argument without his J ≥ N assumption;
- external inputs are Rosser–Schoenfeld's (3.14), θ(x) > x(1 − 1/(2 log x)) for
  x ≥ 563 (needed beyond 2·10^6), genus theory and Hadamard's inequality;
- this is a computer-assisted proof, not yet a formal (Lean) one.

## Next steps

- Adopt Green's deg F₀ cancellation in Lemma 1 (a slightly better criterion)
  and rerun the search.
- An independent human review of NOTES.md before posting.
- Compare Lemma 1 with Green's note, then repeat for the Banks and Kerr
  write-ups and the QRH papers.
- Better optimiser / sharper Archimedean estimates to lower the threshold.

## Running

    cd tenth_disc
    python3 reduction.py                     # table above
    python3 -m unittest test_reduction test_green_explicit test_tail test_rigorous
    python3 certify.py --verify              # re-check the certificate (~1 min)
    python3 rigorous.py                      # interval-arithmetic proof (~5 min)
    python3 certify.py                       # redo the search (~1 hour)

Needs numpy and mpmath; sympy is optional (one symbolic test).

*AI disclosure: the code and write-up were produced with Claude Code (Anthropic).*
