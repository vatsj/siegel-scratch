# Review guide

A checklist for reviewing `NOTES.md`, ordered from "undergraduate-checkable"
to "needs a specialist". For each item, either tick it or write down exactly
where you got stuck. "I could not follow step 3 of Lemma 1" is a useful
review outcome.

Tao's request for AI-generated work: *"do enough due diligence to check the
output (e.g., with an independent AI tool) before posting"*, and disclose the
AI use.

## A. Undergraduate-checkable (linear algebra, elementary number theory)

- [ ] **Lemma 3 (p-adic bound).** It needs Fermat's little theorem, the
  binomial theorem and the product rule.
  - Check (m + n√D)^p ≡ m − n√D (mod p) when D is a non-residue mod p.
  - Check that every derivative of (z1^p − z2)^k Q vanishes mod p^k on E.
  - Check the basis change is unitriangular: compare Green's Lemma 5.1, which
    is word-for-word the same.
- [ ] **Lemma 2 (Hadamard).** Check the Hadamard inequality
  det(Gram) ≤ ∏‖v_i‖². Check the entry bound
  a(a−1)…(a−j+1)·ρ^{a+b−j} ≤ ρ^{a+b}·max(1, a/ρ)^{J−1}. Check
  Σ_{a<R} log⁺(a/ρ) ≤ ∫₀^R log⁺(t/ρ) dt.
- [ ] **M is a positive integer.** Check that Re(u·v̄) = u₁v₁ + d·u₂v₂ for
  u, v ∈ Z[√−d]. Check that a Gram determinant of independent vectors is
  > 0.
- [ ] **Minorant lemma** (in "The inequality (\*)"). Check the formula for
  g(q), that g̃ is nondecreasing, the sum-versus-integral step, and the sign
  of f′.
- [ ] **Bookkeeping of (\*).** Check that the split and ramified corrections
  are *subtracted*. Check that d_lo is used where a bound decreases in d,
  and d_hi where it increases.
- [ ] **Tail algebra.** Plug the constants in and re-derive
  margin ≥ 0.18646 L − 15.609 by hand or in a calculator.
- [ ] **Spot-check by running the code:**
  `python3 -m unittest test_green_explicit` computes exact determinants for
  tiny cases.

## B. Needs first-course algebraic number theory

- [ ] If h(−d) = 1 and d > 163, then d is prime and d ≡ 3 (mod 8). This is
  genus theory, plus the fact that 2 = (a² + 7b²)/4 forces d = 7.
- [ ] A split prime p is the norm of a principal ideal, so
  p = (a² + db²)/4 with b ≠ 0. Hence every p < d/4 is inert.
  - Sanity check: d = 163, where all primes < 41 are inert.
- [ ] The split-prime sum bound in `reduction.split_upper`. This is just
  counting lattice points.

## C. Needs a specialist (or at least careful reading)

- [ ] **Lemma 1 (injectivity), step 3.** The vanishing order of the
  resultant. Compare with Green's Lemma A.1, whose rank proof may be easier
  to follow.
- [ ] **Lemma 1, the rest.** Factorisation over Q and Schwartz–Zippel on the
  rectangle. Compare with Green's Lemma 3.1.
- [x] **Rosser–Schoenfeld.** Theorem 4, eq. (3.14) of *Illinois J. Math.* 6
  (1962), p. 70, states x(1 − 1/(2 log x)) < θ(x) for 563 ≤ x. Checked
  against the paper. Eq. (3.21), used in `reduction.py`, is on the same page.
- [ ] **Code matches maths.** Check that `rigorous.py` computes exactly the
  inequalities stated in NOTES.md.

## Suggested disclosure wording for the blog comment

> AI-generated (Claude Code); all code and notes at <repo link>. The
> argument follows Tao's sketch / Green's note, with every constant made
> explicit. The final inequalities are verified in outward-rounded interval
> arithmetic. Lemmas 1–3 are also checked by exact computation on small
> cases. Human review so far: [describe exactly what was checked, by whom,
> e.g. "sections A–B checked by a reader with an undergraduate maths
> background; Lemma 1 not independently reviewed by a specialist"].
> Independent AI check: [tool, what it was asked, outcome].
