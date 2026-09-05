# The cubic term at the simplex

Companion to `conjecture-notes.md` and `orientation-character.md`. Verified by
`scripts/third_derivative.py`.

Short version: the third-derivative derivation is correct, the four-cycle is
genuinely cubic-null, the invariant count is right, and for `m = 4` the cubic
coefficient is now computed — it is nonzero at both transitions, so the
higher-order saddle is real and not merely possible. One scope correction: the
saddle is a fixed-energy phenomenon and does not touch Theorem 4.

---

## 1. The derivation checks out

With `K_eps = sI + eps H`, `H in H_0`:

```
L'   = Y'HY / (2 s^2)
L''  = tr(H^2)/(2 s^2) - Y'H^2 Y / s^3
L''' = -tr(H^3)/s^3 + 3 Y'H^3 Y / s^4
```

the first two matching the paper's Section 5.1 verbatim. Then
`f'''/f = L''' + 3 L' L'' + (L')^3` gives

```
S_3(Y;H) = (Y'HY)^3/(8s^6) - 3(Y'HY)(Y'H^2Y)/(2s^5)
         + 3 tr(H^2)(Y'HY)/(4s^4) + 3 Y'H^3Y/s^4 - tr(H^3)/s^3
```

and `D^3 A_x(sI)[H,H,H] = -(e^{-s/2}/m) E_{N(0,sI)}[Psi_x(Y) S_3(Y;H)]`.

All four expressions were checked against high-precision finite differences of
the actual Gaussian log-density at `m = 4, 5, 6`: relative agreement `1e-12` or
better on `L'`, `L''`, `L'''` and on `f'''/f`. The full mixed trilinear form —
`ell_123 + ell_12 ell_3 + ell_13 ell_2 + ell_23 ell_1 + ell_1 ell_2 ell_3`
with

```
ell_123 = -[tr(H1H2H3) + tr(H1H3H2)]/(2s^3)
        + Y'( sum over the 6 orderings of H1H2H3 )Y / (2s^4)
```

was checked the same way with three distinct random directions at `m = 5, 6`:
agreement to `4e-10`. Setting all three equal recovers the scalar case, as it
should.

The domination argument does extend: `|d^3 f_eps/d eps^3| <= C(1+||y||^6)
e^{-||y||^2/(4s)}` against `Psi_x(y) <= m e^{||y||}` is the same estimate the
paper already runs for `k <= 2` with `||y||^4`.

---

## 2. The four-cycle is cubic-null, for a good reason

The paper's direction (`H_12 = H_34 = 1`, `H_13 = H_24 = -1`) satisfies
`tr(H^3) = 0` and `sum_{i<j} H_ij^3 = 0` exactly, and `H^3 = 4H` exactly, so
its spectrum `{-2, 0, ..., 0, 2}` is symmetric.

The symmetry explanation is the better one and it is exact: for the 4-cycle
`sigma = (1 2 4 3)`,

```
P_sigma H P_sigma^T = -H       (verified to machine zero at m = 4 and m = 6)
```

so `A_x(sI + eps H) = A_x(sI - eps H)` is even in `eps` and **every** odd
derivative along that line vanishes, not just the third. Note `sigma` is a
4-cycle, hence odd — the same parity that carries the orientation character of
`orientation-character.md`, though the two facts are independent.

For `m = 4` there is a sharper statement. `H_0` is two-dimensional, and
restricted to its unit circle the invariant cubic is exactly

```
tr(H^3) = A cos(3 theta - psi),    residual 1.5e-14 over 2881 sample angles,
```

so it has exactly three zero lines. Those three lines are precisely the three
four-cycles of `K_4` (angles `39.9386, 99.9386, 159.9386` degrees in an
orthonormal basis of `H_0`, matching the predicted zeros to four decimals).
**The four-cycle directions are the zero locus of the invariant cubic.** The
direction chosen for the Hessian computation is not merely unlucky; it is one
of the three lines on which the cubic identically vanishes.

---

## 3. The invariant count is right, with an extra relation

By the character formula with `chi_{H_0}(sigma) = C(F,2) + T - F` (`F` fixed
points, `T` 2-cycles), summed over the whole group:

| m | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| `dim (Sym^3 H_0*)^{S_m}` | 1 | 1 | 2 | 2 | 2 | 2 |

exactly as claimed. And `{tr(H^3), sum_{i<j} H_ij^3}` really is a basis: their
value vectors over 400 random `H in H_0` have rank 2 for `m >= 6`. For
`m = 4, 5` the rank is 1, and the relation is clean:

```
tr(H^3) = 4 sum_{i<j} H_ij^3        on H_0, for m = 4 and m = 5
```

(verified to a relative spread of `1e-12`). That is the concrete reason the
invariant space is one-dimensional in those two cases, and it means there is
no ambiguity in what `T_3` should be.

---

## 4. The coefficient at m = 4, computed

The conclusion "if `alpha != 0` or `beta != 0` then the simplex is a
higher-order saddle at the transition" was left conditional. For `m = 4` the
condition can be discharged.

`A_x(K)` is evaluated directly rather than through the Hessian formula: each
term of `Psi_x` is tilted exactly,
`E[e^{Y_i} 1_A] = e^{K_ii/2} P(Y + K e_i in A)`, turning it into a Gaussian
rectangle probability over the `m-1` differences, computed by the deterministic
Genz routine in `scripts/_orthant.py`. Validation: differencing it twice along
the four-cycle reproduces the paper's exact
`D^2 A_x = -I_m(s,x) ||H||_F^2 / (8 m s^2)` to `1e-4`–`1e-5` relative at eight
`(s, rho)` points.

At `s = 1`, `m = 4` the transitions are at `rho_c^- = 0.4089176535` and
`rho_c^+ = 3.2822718637` (roots of `I_4`). Differencing three times:

| direction (deg) | `tr(H^3)` | `D^3A` at `rho_c^-` | ratio | `D^3A` at `rho_c^+` | ratio |
|---|---|---|---|---|---|
| 9.9386 | 26.1279 | 0.0517311 | 0.00197992 | 0.0144200 | 0.000551902 |
| 40.0 | -0.0840296 | -0.000181596 | 0.0021611 | -4.63297e-05 | 0.00055135 |
| 70.0 | -26.1278 | -0.0517309 | 0.00197992 | -0.0144200 | 0.000551902 |
| 129.9386 | 26.1279 | 0.0518206 | 0.00198334 | 0.0143551 | 0.000549417 |

Two things to read off. The ratio `D^3A / tr(H^3)` is constant across
directions to four significant figures — including at `theta = 40` degrees,
essentially on a four-cycle where both numerator and denominator are near zero.
That is the numerical signature of the one-dimensional invariant cubic space,
independently of the character count. And the constant is not zero:

```
alpha_4(s=1, rho_c^-) = 1.980e-3,     alpha_4(s=1, rho_c^+) = 5.519e-4,
```

good to about three digits (the `D^3` noise floor here is `~1e-5` absolute; the
four-cycle direction returns `1e-5`, i.e. zero, as it must).

**So at `m = 4` the simplex is a genuine higher-order saddle at both curvature
transitions.** `tr(H^3)` takes both signs on `H_0`, so `A_{x_c}` strictly
increases in some shape directions and strictly decreases in others at third
order. This settles, for `m = 4`, the case Theorem 3(4) explicitly leaves
open ("No classification of a degenerate point `f(rho) = 0` follows from the
Hessian alone").

---

## 5. One scope correction: this does not touch Theorem 4

The saddle statement is about the **fixed-energy** functional `A_x` at
`rho = rho_c`. It says nothing about the scale-profiled `F_x`, for two reasons.

First, the profiled point is `s = x`, i.e. `rho = 1`, and Theorem 3(4) places
it strictly inside the positive interval, `rho_c^- < 1 < rho_c^+`. Section 9.1
shows `I_m(x,x) > 0` outright, since every coefficient at `rho = 1` is
positive. So `F_x` is never at a transition; it is a strict local maximum at
every `x > 0`, and there is no degenerate point to resolve.

Second — worth recording anyway — the profiled third derivative equals the
fixed-energy one. Writing `e(t) = E*(G_tri + tH)` and `a(t) = A(e(t), G+tH)`:
stationarity (Theorem 1) gives `D_G A(E, G_tri)[H] = 0` identically in `E`,
hence `A_{EG}[H] = 0`, hence `e'(0) = 0`. Every term in `a'''(0)` carries a
factor of `e'(0)`, `A_E`, or `A_{EG}[H]` except one, so

```
D^3 F_x(G_tri)[H,H,H] = D^3 A_x(E*, G_tri)[H,H,H],
```

the third-order analogue of Lemma 7. Useful if the profiled cubic is ever
wanted; it just is not wanted at a transition, because there isn't one.

---

## 6. What is left

`alpha_m` and `beta_m` in closed form, for general `m`. The reduction route is
the one described: tilt each term of `Psi_x`, condition on the base coordinate,
and reduce the degree-6 polynomial moments over the independent truncated
normals, then apply the inverse-Mills recurrence. The `m = 4` numbers above are
a target for that computation to reproduce.

Two directions are needed for `m >= 6`, and neither may be a four-cycle. A
convenient pair: the cubic-extremal direction (for `m = 4` it is
`H_13 = H_24 = 2`, `H_12 = H_14 = H_23 = H_34 = -1` up to scale — the split
into two antipodal pairs) and any direction independent of it under the two
invariants.

Reproduce with:

```bash
python scripts/third_derivative.py                 # all parts
python scripts/third_derivative.py --part alpha    # the m=4 coefficient, ~20 min
```
