# Green's argument, made explicit, applied to the tenth discriminant

*Working notes for Tao's
[AI-math challenge](https://terrytao.wordpress.com/2026/10/20/an-experimental-ai-math-challenge-using-tenth-discriminant-bounds-to-explore-arguments-ruling-out-siegel-zeroes/).
They follow Tao's [sketch of Green's argument](https://terrytao.wordpress.com/2026/10/19/a-sketch-of-the-proof-of-non-existence-of-siegel-zeroes/).
Green's own note, and the Banks/Kerr preprints, were not reachable from the
environment these notes were written in. So the general injectivity lemma
below is my own completion of the "irreducible P" cheat in the sketch, not
Green's argument.*

**AI disclosure:** produced with Claude Code (Anthropic). Every lemma is
checked by brute force on small cases (`test_green_explicit.py`). The
certificate is re-checked in outward-rounded interval arithmetic
(`rigorous.py`, mpmath.iv at 160 bits).

## Result

> **Claim (verified in interval arithmetic).** Green's polynomial-method argument, with
> the explicit constants below, shows that an imaginary quadratic field
> Q(√−d) with d > 163 and class number one must have d < 515,960.
> An elementary check (every prime d ≡ 3 mod 8 in (163, 515,960] has a prime
> p < d/4 with (−d/p) = +1) disposes of the rest. So the method recovers
> the class number one theorem without Heegner/Stark/Baker, modulo the
> caveats at the end.

For comparison: Goldfeld–Oesterlé give d ≤ e^55 ≈ 7.7·10^23.

| range of d | how it is excluded | where |
|---|---|---|
| (163, 515,960] | some p < d/4 splits, so h > 1 (elementary) | `certify.small_d_check` |
| [515,960, 10^200] | finitely many intervals, each with explicit parameters (N1, N2, R, S, J) violating (\*) | `results/certificate.json`, `certify.py --verify` |
| [10^40, ∞) | closed form: margin ≥ 0.18 log d − 16 > 0 | `tail.py`, derivation below |

## Setup

Assume h(−d) = 1 with d > 163. By genus theory d is prime and d ≡ 3 (mod 8).
A split prime is the norm of a principal ideal, p = (a² + d b²)/4 with b ≠ 0.
Hence:

* every prime p < d/4 is inert (except that p = d ramifies);
* split primes below x satisfy Σ log p/p ≤ 4 log x · ((π²/6)/d + (π/√d)(1 + log(2√(x/d))))
  (`reduction.split_upper`).

Put D = −d. Parameters are a grid N1 × N2, degrees (R, S) and a derivative
order J.

* Λ = { P ∈ Z[z1, z2] : deg_{z1} P < R, deg_{z2} P < S }, a lattice of rank RS.
* E = { (α, ᾱ) : α = m + n√D }, with m running over N1 consecutive integers centred
  at 0 and n over N2. So |α|² ≤ ρ² := (N1//2)² + d (N2//2)².
* T P = (∂_{z1}^j P(e))_{j<J, e∈E} ∈ Z[√D]^{J·|E|}.
* M = Gram determinant of T Λ, using the real inner product Re(u v̄).
  For u, v ∈ Z[√D] this is u₁v₁ + d u₂v₂ ∈ Z, so M ∈ Z.

## Lemma 1 (injectivity), explicit and general

**If J (N1 N2 − S + 1) > 2(R−1)(S−1) + max(N1, N2)(R + S − 2), then T is
injective on Λ, so M is a positive integer.**

*Proof.* Let P ≠ 0 vanish to horizontal order ≥ J at every point of E. Factor
P = A(z2) ∏ F_i^{e_i}, where the F_i are distinct irreducibles of z1-degree
r_i ≥ 1 and z2-degree s_i. Then Σ e_i r_i ≤ R − 1 and deg A + Σ e_i s_i ≤ S − 1.

1. The z2-coordinates ᾱ of E are pairwise distinct, and A ≠ 0 has degree ≤ S−1.
   So A vanishes at ≤ S − 1 points of E. Let E' be the remaining points, so
   |E'| ≥ N1N2 − S + 1.
2. Take e = (z0, c) ∈ E'. Then F_i(·, c) ≢ 0: otherwise (z2 − c) and its
   conjugate would divide F_i, contradicting irreducibility with r_i ≥ 1.
   Let μ_i(e) be the order of z0 as a root of F_i(·, c). The vanishing
   hypothesis gives Σ_i e_i μ_i(e) ≥ J.
3. Res_{z1}(F_i, ∂F_i) is nonzero (F_i is squarefree in Q(z2)[z1], char 0)
   and has z2-degree ≤ (2r_i − 1)s_i. Expand the Sylvester determinant
   multilinearly, writing each row as (z1 − z0)^{μ−1}(…) + (z2 − c)(…), as in
   Tao's sketch. This shows (z2 − c)^{μ_i(e) − 1} divides the resultant. The
   c's are distinct, so Σ_e (μ_i(e) − 1)_+ ≤ (2r_i − 1)s_i ≤ 2 r_i s_i.
4. G(m, n) := F_i(m + n√D, m − n√D) is a nonzero polynomial: the linear change
   of variables is invertible. Its total degree is ≤ r_i + s_i. By
   Schwartz–Zippel on the N1 × N2 box, F_i vanishes at
   ≤ (r_i + s_i) max(N1, N2) points of E.

Now combine. Σ_{e∈E'} Σ_i e_i μ_i(e) ≤ Σ_i e_i [(μ−1)_+ count + zero count].
This is ≤ 2(R−1)(S−1) + max(N1,N2)(R + S − 2), using
Σ e_i r_i s_i ≤ (R−1)(S−1). The left side is ≥ J|E'|, a contradiction. ∎

## Lemma 2 (Archimedean)

Use Hadamard's inequality with the monomial basis z1^a z2^b. Every entry of
T(z1^a z2^b) is a^{(j)} α^{a−j} ᾱ^b, and a^{(j)}/ρ^j ≤ max(1, a/ρ)^{J−1}.
Since Σ_{a<R} log⁺(a/ρ) ≤ R log(R/ρ) − R + ρ, this gives

    log M ≤ RS log(|E| J) + RS(R + S − 2) log ρ + 2S(J − 1)(R log(R/ρ) − R + ρ).

## Lemma 3 (non-Archimedean)

Let p be an odd prime with (D/p) = −1. Then (m + n√D)^p ≡ m − n√D (mod p), so
z1^p − z2 ∈ pZ[√D] on E. Every z1-derivative of z1^p − z2 is divisible by p,
so every z1-derivative of (z1^p − z2)^k Q is ≡ 0 (mod p^k) on E. Consider the
basis z1^{i−pk}(z1^p − z2)^k z2^j with k = k_{ij} = min(⌊i/p⌋, S−1−j). It is
unitriangular with respect to z1^i z2^j, so

    v_p(M) ≥ 2 K_p,    K_p = Σ_{i<R, j<S} min(⌊i/p⌋, S−1−j)   (closed form in `K_exact`).

Brute force shows this bound is often sharp (v_p(M) = 2K_p exactly).

## The inequality (\*) and how it is evaluated

If h(−d) = 1, then

    2 Σ_{p odd inert} K_p log p  ≤  log M  ≤  (Lemma 2).          (*)

The left side is bounded below as follows:

* Take the sum over **all** odd p < R. Subtract the split primes using
  K_p ≤ S R(R−1)/(2p) together with the split-prime bound, and subtract the
  ramified prime p = d if d < R.
* For p ≤ 2·10^6, sum K_p exactly.
* Beyond that, use Abel summation against the decreasing minorant
  f(t) = t·G((R−1)/t − 1) ≤ K_t, where G is the integral of the continuous
  extension of q ↦ Σ_j min(q, S−1−j). Combine it with Rosser–Schoenfeld,
  θ(t) > t(1 − 1/(2 log t)) for t ≥ 563, on a fine geometric grid. Only
  monotonicity is used, so the discretisation is one-sided.

If some admissible parameters make (\*) fail for all d ∈ [d_lo, d_hi], none of
those d has class number one. The split and ramified terms are decreasing in
d, so they are evaluated at d_lo. ρ is increasing in d, so it is evaluated at
d_hi.

## The tail d ≥ 10^40 (closed form)

Take R = ⌊d/4⌋, S = ⌊R/64⌋, J = 2S, N2 = ⌈d^{1/4}/2⌉ and N1 = ⌈1.1R/N2⌉.
Every odd p < R is then < d/4, so it is inert, and no split or ramified
primes occur.

* **Injectivity.** S ≥ R/65 gives (R + S)/(2S) ≤ 33.1. So Lemma 1 holds as
  soon as N1(N2 − 33.1) > 1.016R, which is true since N1N2 ≥ 1.1R and
  N2 > 10^9.
* **Prime side.** Drop p ≤ T = 10^6 and set η = 1/(2 log T). Then

      Σ_{p<R} f(p) log p ≥ (1−η)(∫₀^{R−1} f − ∫₀^T f) − f(T)(θ(T) − (1−η)T).

  Here ∫₀^{R−1} f = (R−1)² S·I(S) with the **exact** identity
  I(S) = ½ log S − ½ + (log S + 2)/(4S) (checked with sympy). Also
  ∫₀^T f ≤ (R²S/2)(1 + log(TS/R)), f(T) ≤ R²S/(2T), and θ(10^6) = 998484.175….
  With L = log d and log S ≥ L − log 256 − o(1), the normalised prime side is
  ≥ 0.96381 L − 16.614.
* **Archimedean side.** ρ² ≤ 0.138125 d^{3/2}(1 + 10^{−9}), so
  log ρ ≤ ¾L − 0.98977. Also log R − log ρ ≤ L/4. The normalised Lemma 2
  bound is then ≤ (65/64)(¾L − 0.98977) + L/64 + o(1) = 0.77734 L − 1.0052.

So the margin is ≥ 0.18646 L − 15.609 − o(1) ≥ 0.18 L − 16 > 0 for L ≥ log 10^40.
`rigorous.check_tail` recomputes these constants in interval arithmetic and
absorbs every o(1) term into 10^{−38} L + 10^{−8}. Each such term is
log(1 − 4/d), log(1 − 64/R), (R−1)²/R², 3L/R, or a 10^{−9} correction in the
bounds for log ρ.
`tail.py` evaluates the same closed form without the linearisation. The
tests check it against the full numerical engine.

## Caveats

* **Rigour.** `rigorous.py` re-evaluates every certificate interval and the
  tail constants with outward-rounded interval arithmetic. Floats only choose
  the Abel grid points. All 18 margins are certified positive (smallest
  0.00172, in units of R²S), and they agree with the float values to 6
  digits. The script takes about 5 minutes.
* Lemma 1 in the general (reducible P) case is my own argument. It should be
  compared against Green's note.
* External inputs: Rosser–Schoenfeld (1962, Thm 4),
  θ(x) > x(1 − 1/(2 log x)) for x ≥ 563, used only for x ≥ 2·10^6 and checked
  directly below that (`test_rigorous.py`); genus theory; the norm form of
  split primes; Hadamard's inequality.
* The parameter search is heuristic, so a better optimiser may lower 515,960.
