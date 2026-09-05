# The quartic term along the four-cycle

Companion to `third-derivative.md`. Verified by `scripts/quartic_four_cycle.py`.

Short version: the derivation is correct in every step I can check, including
the closed forms for `C` and `D` and the final constant, and the numerical
table reproduces digit for digit. The one thing to correct is the
interpretation: the four-cycle slice is exactly the slice on which the cubic
vanishes, so the quartic normal form and the square-root law describe that
slice and not the local geometry at the transition, which is a cubic saddle.

---

## 1. Verified

**The reduction to one correlation.** `H = u v' + v u'` with `u = e_1 - e_4`,
`v = e_2 - e_3`; `u'v = 0`, `|u|^2 = |v|^2 = 2`, `u'Hu = v'Hv = 0`,
`u'Hv = 4`, and `u'Hw = 0` for `w` orthogonal to both. So along `K_eps`, the
standardized pair `(U,V)` has correlation `r = 2 eps/s` and nothing else moves;
`(U,V)` is independent of the orthogonal block. That is genuinely why the
four-cycle factorizes.

**Mehler and the fourth score.** `d^k/d eps^k (f_eps/f_0)|_0 =
(2/s)^k He_k(U) He_k(V)`, giving `(P-2s)(Q-2s)/s^4` at `k = 2` — the paper's
own Hessian integrand — and

```
(P^2 - 12sP + 12s^2)(Q^2 - 12sQ + 12s^2) / s^8
```

at `k = 4`. Checked against high-precision finite differences of the true
Gaussian density at `m = 6`, `s = 1.3`: `7e-10` at `k = 2` and `3e-8` at
`k = 4`. The general route agrees: `L'''' = 3 tr(H^4)/s^4 - 12 Y'H^4Y/s^5`
matches finite differences to `6e-9`, and for this `H`, `tr H^2 = 8`,
`tr H^3 = 0`, `tr H^4 = 32`, `H^3 = 4H`, `H^4 = 4H^2` all hold exactly. The
Bell-polynomial assembly `f''''/f = L'''' + 4L'L''' + 3(L'')^2 + 6(L')^2L'' +
(L')^4` is the right expansion, and the corrected `6(L')^2 L''` term is what
makes it collapse to the Hermite product.

**The truncated moments.** I re-derived `C` and `D` by hand from
`mu_1 = -lam`, `mu_2 = 1 - t lam`, `mu_3 = -(t^2+2) lam`,
`mu_4 = 3 - (t^3+3t) lam`. Both reproduce exactly:

```
E Z^2 = 2 - 2t lam - 2 lam^2
E Z^4 = 12 - 2t^3 lam - 18t lam - 2t^2 lam^2 - 16 lam^2
D     = E Z^4 - 12 E Z^2 + 12 = -2 lam [ t^3 - 3t + lam(t^2 - 4) ]
C     = He_4(c) + lam [ 4c^3 - 6c^2 t + 4c t^2 - 16c - t^3 + 9t ]
```

**The counting and the constant.** `D(t)` is `4 E[He_4((X-X')/sqrt 2)]` and
`C(t)` is `4 E[He_4((c-X)/sqrt 2)]`; the two factors of 4 cancel the `16/s^4`
from `(2/s)^4`, which is why the final constant comes out as the clean
`-K_m/(m s^4)`. Four bases lie in `{1,2,3,4}` and give `C*D`, the rest `D^2`.
`He_4` is even, so it does not matter which end of a pair carries the base.

**End to end.** `-K_m/(m s^4)` against direct fourth differences of `A_x`
(evaluated by tilting each `Psi_x` term into a Gaussian rectangle probability,
the same evaluator validated against the exact Hessian to `1e-4`):

| m | s | rho | `-K_m/(m s^4)` | 4th difference, `h -> 0` | rel |
|---|---|---|---|---|---|
| 4 | 1 | 0.6 | -1.74173682 | -1.74312706 | 8.0e-4 |
| 4 | 1 | 1.5 | -0.844327324 | -0.844187255 | 1.7e-4 |
| 5 | 1 | 1.0 | -1.25913721 | -1.25964152 | 4.0e-4 |
| 4 | 2 | 1.2 | -0.0195595243 | -0.0195542864 | 2.7e-4 |
| 6 | 1 | 0.8 | -0.945068041 | -0.949698898 | 4.9e-3 |

One methodological note, since it cost me an hour: the fourth difference must
be extrapolated on *small* steps only. Fitting `h^2 + h^4` across
`h = 0.20 ... 0.08` gives `-1.756` and a spurious 0.8% discrepancy; the `h^6`
term is not negligible at `h = 0.2`. Restricted to the three smallest steps the
fit lands on `-1.74313`, and the residual is stable under refining the orthant
integrator, so it is extrapolation error, not a defect in `K_m`.

**Section 11's table.** Reproduced exactly by an independent implementation of
`C`, `D`, `W_b` and the counting — all ten numbers, `s = 1` and `s = 4`,
`m = 4, 6, 10`, including `rho_c^- = 0.408918`, `rho_c^+ = 3.282272` at
`m = 4, s = 1`, which also match the roots found independently in
`third-derivative.md`.

**Evenness.** `sigma = (2 3)` gives `P_sigma H P_sigma' = -H` (verified to
machine zero) — simpler than the 4-cycle `(1 2 4 3)` used earlier, and the same
conclusion: `eps -> A_x(sI + eps H)` is exactly even, so all odd derivatives
vanish along this line.

---

## 2. The correction: this is the cubic-null slice

Sections 10 and 12-14 read the quartic as *the* local classification at the
transition, and derive `eps ~ |rho - rho_c|^{1/2}` as *the* throat law. Both
are true on the four-cycle slice and false off it.

From `third-derivative.md`: at the transition the Hessian vanishes on all of
`H_0`, but the cubic does not. For `m = 4`, `s = 1` the coefficient is
`alpha_4 = 1.980e-3` at `rho_c^-` and `5.519e-4` at `rho_c^+`, both nonzero.
So the simplex at a transition is a **cubic saddle** in `H_0`, and the leading
correction in a generic direction is third order, not fourth.

The four-cycle is precisely the exception. It annihilates both invariant
cubics, and for `m = 4` the three four-cycles of `K_4` are *exactly* the zero
lines of the invariant cubic on the two-dimensional `H_0`. The evenness
argument of Section 9 is the reason: it forces the odd part to vanish on that
line, which is the same thing as saying the line lies in the cubic's zero
locus. So the fourth-order analysis is not seeing the generic behaviour at the
transition; it is seeing what is left after the symmetry of this particular
direction has removed the leading term.

The bifurcation scalings differ accordingly. Write `mu = rho - rho_c` and
`a = -I'_rho |H|_F^2/(8 m s^2)`.

* On the four-cycle, `A - A_* = a mu eps^2 + b eps^4` with
  `b = -K_m/(24 m s^4)`, stationary branches at `eps^2 = -(a/2b) mu`, so
  **`eps ~ |mu|^{1/2}`** — the square-root law, correct as derived.
* In a generic direction `Hhat` with `T_3(Hhat) != 0`,
  `A - A_* = a mu eps^2 + c eps^3` with `c = alpha_m tr(Hhat^3)/6`, stationary
  at `eps = -2a mu/(3c)`, so **`eps ~ |mu|`**.

Since `|mu| << |mu|^{1/2}` for small `mu`, the generic branch sits *closer* to
the simplex than the four-cycle branch: the throat is narrowest in the
directions the four-cycle cannot see. Both scalings do regularize the divergent
quadratic horn — the level set `Delta A = const` gives `eps ~ (Delta A/b)^{1/4}`
on the slice and `eps ~ (Delta A/c)^{1/3}` generically — so the conclusion
"not Gabriel's horn" stands. It is the exponent that is slice-dependent, and
the honest statement is a *cubic* throat generically, with a quartic throat on
the cubic-null locus.

---

## 3. The profiled fourth derivative, explicitly

Section 15 is right that the scale correction no longer drops out at fourth
order. It can be written down. With `e(t) = E*(G + tH)`, stationarity gives
`A_E = 0`, `A_EG[H] = 0` (the paper's Theorem 1 makes `D_G A(E, G_tri)[H]`
vanish identically in `E`), hence `e'(0) = 0` and
`e''(0) = -A_EGG[H,H]/A_EE`. Collecting the surviving terms,

```
D^4 F_x(G_tri)[H^4] = D^4 A_x[H^4] + 6 A_EGG e''(0) + 3 A_EE e''(0)^2
                    = D^4 A_x[H^4] - 3 (A_EGG[H,H])^2 / A_EE .
```

`A_EE < 0` at the profiling maximum, so the correction is positive: the
profiled quartic always exceeds the fixed-energy one. (At third order the same
collection leaves `D^3 F = D^3 A`, as recorded in `third-derivative.md`.)

This is moot for the transition, though, and for the reason given there: the
profiled point is `s = x`, i.e. `rho = 1`, which Theorem 3(4) places strictly
inside `(rho_c^-, rho_c^+)`, and Section 9.1 shows `I_m(x,x) > 0` outright. The
profiled functional is never at a transition.

---

## 4. The conjectured sign law is false

`K_m(s, rho_c^- s) > 0` held in every case I tested. **`K_m(s, rho_c^+ s) < 0`
does not.** It fails for all `m` past a finite threshold:

| s | 0.25 | 1 | 4 | 9 | 25 |
|---|---|---|---|---|---|
| `K_m(rho_c^+) < 0` holds up to | m = 11 | m = 10 | m = 11 | m = 25 | m = 59 |
| first `m` with `K_m(rho_c^+) > 0` | 12 | 11 | 12 | 26 | 60 |

At `s = 1`: `K(9) = -0.02741`, `K(10) = -0.007934`, `K(11) = +0.007972`,
`K(12) = +0.02112`, rising to `+0.0492` at `m = 15`, `+0.0878` at `m = 25`,
`+0.0984` at `m = 60`. Two independent runs give the same threshold. This is
not a quadrature artifact: a peak-centred quadrature agrees with the
`b`-centred one to every printed digit at `m = 10 ... 40`, and the roots are
genuine sign changes of `I_m` (`+0.001048` just below, `-0.0008679` just above,
at `m = 25, s = 1`, with `I_m` still negative out at `rho = 10`).

**There is a structural reason, and it does not depend on the numerics.** Split

```
K_m = (m-4) * Idd  +  4 * Icd,     Idd = int W_b D^2,   Icd = int W_b C D.
```

`Idd > 0` always: the integrand is a positive weight times a square, and `D` is
not identically zero. So `K_m` carries a strictly positive part that grows
linearly in `m`, and the law can only survive while `Icd` is negative enough to
cancel it. At `rho_c^+` it is not, for two compounding reasons — `Idd` shrinks
only slowly while `(m-4) Idd` grows, and `Icd` itself changes sign:

| m (s = 1) | `rho_c^+` | `Idd` | `Icd` | `(m-4) Idd` | `K_m` |
|---|---|---|---|---|---|
| 6 | 3.404000 | 0.00663072 | -0.0338309 | 0.0132614 | -0.122062 |
| 10 | 3.571773 | 0.00433014 | -0.00847865 | 0.0259809 | -0.00793374 |
| 12 | 3.634812 | 0.0037215 | -0.00216204 | 0.029772 | +0.0211238 |
| 15 | 3.713795 | 0.00308713 | **+0.00382351** | 0.0339585 | +0.0492525 |
| 25 | 3.900524 | 0.00198007 | +0.0115577 | 0.0415815 | +0.0878123 |
| 60 | 4.231599 | 0.000838813 | +0.0128544 | 0.0469735 | +0.0983911 |

Past `m = 12` or so both terms are positive and there is nothing left to make
`K_m` negative. So no reduction through the Mills recurrence is going to prove
`K_m(rho_c^+) < 0`; the statement is simply not true. The `rho_c^-` half looks
robust for the same reason read forwards — there `Icd > 0` as well, so both
terms push the same way.

Two consequences for Sections 10 and 13.

* For `m >= 11` (at `s = 1`), `K_m > 0` at *both* transitions, so `D^4 A < 0`
  at both: the simplex is quartically **maximizing** along the four-cycle at
  `rho_c^+` too, not minimizing. The "mirror versions in curvature" reading of
  the two crossings holds only for small `m`.
* The bifurcation side flips with it. At `rho_c^+` we have `a > 0` (since
  `I'_rho < 0` there); once `K_m > 0`, `b = -K_m/(24 m s^4) < 0`, so `a/b < 0`
  and `eps^2 = -(a/2b) mu > 0` requires `mu > 0` — the opposite side from the
  small-`m` case.

What survives is the derivation, which is exact, and the small-`m` picture. The
right open question is not "prove the sign law" but "for which `(m,s)` is
`K_m(rho_c^+)` negative" — the threshold curve `m*(s)`. It is not monotone at
the small-`s` end (12, 11, 12 at `s = 0.25, 1, 4`) but grows steadily after
that (26 at `s = 9`, 60 at `s = 25`), so the law is worst at moderate `s` and
the four-cycle quartic keeps its small-`m` sign longer as the drift grows. The
`s = 25` row sits at `|K| ~ 1e-17` in absolute terms; the crossing there is
clean in the quadrature but is the least comfortable row in the table.

## 5. Also open

The generic quartic (not just along the four-cycle) would need the full
`Sym^4` invariant decomposition of `H_0`, which has more than the two
invariants that appear at third order.

Reproduce with:

```bash
python scripts/quartic_four_cycle.py                    # all parts
python scripts/quartic_four_cycle.py --part finite-diff # the end-to-end check
python scripts/quartic_four_cycle.py --part signs       # the failure and its cause
```
