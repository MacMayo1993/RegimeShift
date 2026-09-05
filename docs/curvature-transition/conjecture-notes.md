# Notes on the two open conjectures

Working notes on Conjecture 15 (exactly two curvature zeros) and Conjecture 16
(global scale-profiled cyclic maximality) of *A curvature transition at the
regular simplex*, rev. 6. Everything numerical here is reproduced by
`scripts/curvature_conjectures.py`; everything proved is proved below.

Notation follows the paper: `n = m - 1`, `b = sqrt(s)(1 + rho)`,
`R = phi/Phi`, and

```
J_j(b) = integral phi(t-b) phi(t)^j Phi(t)^{n-j} dt,   j = 1,2,3.
```

---

## 0. The normalized curvature is the right object

Dividing (1.11) by the strictly positive `4 sqrt(s) J_1(b)` gives

```
Q_{m,s}(rho) = (3 rho - 1)/2  -  (s/4)(1 + rho)(1 - rho)^2
             + (sqrt(s)/12)[ n(3 rho^2 + 10 rho - 5) - 3(1 - rho)^2 ] u(b)
             + (n(n+1)/6)(3 rho - 1) v(b),
             u = J_2/J_1,  v = J_3/J_1,
```

which has exactly the zeros and signs of `I_m(s, rho s)`. This matches
Theorem 2 term by term. As a check on the whole chain, the closed
three-integral form (1.11) was compared against the raw Hessian integrand
(1.8) at 40-digit precision: they agree to at least 33 digits at every test
point (`--part theorem2`).

---

## 1. Conjecture 15: the concavity route survives an adversarial search

Write `Q = A + B u + C v` with

```
A = (3 rho - 1)/2 - (s/4)(1 + rho)(1 - rho)^2
B = (sqrt(s)/12)[ n(3 rho^2 + 10 rho - 5) - 3(1 - rho)^2 ]
C = (n(n+1)/6)(3 rho - 1)
```

Since `db/drho = sqrt(s)`,

```
Q'' = A''  +  B'' u  +  2 sqrt(s) B' u'  +  s B u''
           +  2 sqrt(s) C' v'  +  s C v''                        (*)
```

with `' = d/db` on `u, v` and `' = d/drho` on `A, B, C`. The coefficient
derivatives are elementary and worth recording, because two of them settle
signs outright:

```
A''  = -(s/2)(3 rho - 1)
B'   = (sqrt(s)/12)[ n(6 rho + 10) + 6(1 - rho) ]
B''  = sqrt(s)(n - 1)/2
C'   = n(n+1)/2,        C'' = 0
```

**The polynomial part of `Q''` is exactly `-(s/2)(3 rho - 1)`.** It is
negative precisely on `rho > 1/3` and vanishes at `rho = 1/3`. This is a clean
structural reason the target is stated on `rho > 1/3`, and it also says that at
`rho = 1/3` the sign of `Q''` is carried entirely by the integral terms, with
no polynomial help. That is the sharpest point of the conjecture, so the scan
below samples `rho` within `1e-9` of `1/3`.

Signs in (*), given `u, v > 0` and Proposition 14 (`u' , v' < 0`):

| term | sign | why |
|---|---|---|
| `A'' = -(s/2)(3 rho - 1)` | **negative** on `rho > 1/3` | polynomial |
| `B'' u = sqrt(s)(n-1) u / 2` | **positive** | `n >= 3`, `u > 0` |
| `2 sqrt(s) B' u'` | negative | `B' > 0` for `rho >= 1/3, n >= 3`, since `n(6 rho + 10) + 6 - 6 rho >= 12 rho + 36 > 0`; and `u' < 0` |
| `s B u''` | unknown | `u''` not signed |
| `2 sqrt(s) C' v' = n(n+1) sqrt(s) v'` | negative | `C' > 0`, `v' < 0` |
| `s C v''` | unknown | `C >= 0` on `rho >= 1/3`; `v''` not signed |

So there is exactly **one** term that is unconditionally positive, `B'' u`,
plus two of unknown sign. Any proof has to dominate `sqrt(s)(n-1) u / 2`.

### Where the domination comes from

Fix `rho` and let `s -> 0`, so `b -> 0`. Then `A''` is `O(s)` while both
`B'' u` and `2 sqrt(s) C' v'` are `O(sqrt(s))`. The polynomial term is
asymptotically irrelevant and

```
Q''(rho) / sqrt(s)  ->  L(n) = (n-1) u(0)/2  +  n(n+1) v'(0),
```

uniformly in `rho`. The sign of `Q''` at small `s` is therefore decided by a
single constant per `n`, and it is the `n(n+1) v'(0)` term — not the
polynomial — that beats the positive term:

| n | u(0) | v'(0) | (n-1)u(0)/2 | L(n) |
|---|---|---|---|---|
| 3 | 0.535596 | -0.198219 | 0.535596 | **-1.84303** |
| 4 | 0.458610 | -0.139085 | 0.687915 | **-2.09378** |
| 5 | 0.403146 | -0.103873 | 0.806293 | **-2.30989** |
| 10 | 0.259145 | -0.0387164 | 1.16615 | **-3.09266** |
| 30 | 0.116461 | -0.00675634 | 1.68869 | **-4.59471** |
| 100 | 0.0441657 | -0.000850423 | 2.18620 | **-6.40307** |
| 1000 | 0.00594859 | -1.26759e-5 | 2.97132 | **-9.71727** |
| 10000 | 0.000723278 | -1.61734e-7 | 3.61603 | **-12.5589** |

`L(n)` is comfortably negative and appears to grow like `-c log n`, while the
positive term grows only like `(n-1)u(0)/2 ~ sqrt(2 log n)/2`. This is the
first sub-inequality a proof should target:

> **Target 1.** `(n-1) u(0)/2 + n(n+1) v'(0) < 0` for every `n >= 3`.

It is a statement about `b = 0` only — two explicit integrals — and for `n = 3`
the paper's elementary formulas for `J_2, J_3` apply directly.

### An exact identity that makes `Q''` computable in closed moment form

Under `p_b(t) = phi(t-b) phi(t) Phi(t)^{n-1} / J_1(b)`,
`(log p_b)'(t) = b - 2t + (n-1) R(t)`. Integrating by parts (boundary terms
vanish by Gaussian decay, `R` grows only linearly at `-infinity`) gives, for
suitable `f`,

```
E_b[f'(T)] = -b E_b[f] + 2 E_b[T f] - (n-1) E_b[f R].
```

Taking `f == 1` yields

> **(A)**  `E_b[T] = b/2 + (n-1) u(b)/2`,

and substituting (A) back gives

> **(B)**  `Cov_b(T, f(T)) = E_b[f'(T)]/2 + ((n-1)/2) Cov_b(f(T), R(T))`.

Both were verified to 26+ significant digits at `n = 3, 5, 9, 100` and
`b = -1.1, 0, 1.3, 2.7, 4` (`--part identities`).

(B) is more useful than plain `d/db E_b f = Cov_b(T, f)` because it trades the
`T`-covariance for an `R`-covariance plus a derivative, and `R' = -R(T+R)` is
explicit. For instance, with `f = R^2`,

```
v' = -E_b[R^2 (T + R)]  +  ((n-1)/2) Cov_b(R^2, R).
```

Here `T + R(T) = E[T - X | X <= T] > 0` (the paper's own observation), so the
first term is negative, while `Cov_b(R^2, R) > 0` because `R` is strictly
decreasing and hence `R` and `R^2` are comonotone. So (B) exhibits `v'` as a
difference of two explicitly signed quantities — exactly the "covariance
representation whose combined sign is controlled" the proof plan asks for. It
also re-derives (A) as the `f == 1` case, which is the identity that makes the
small-`s` limit above computable in closed form.

### The scan

`Q''` was evaluated at 30-digit precision on a deliberately hostile grid:
`n` in {3, 4, 5, 8, 20, 100, 1000, 9999}, `s` in {1e-6, 1e-3, 0.03, 0.3, 1, 3,
10, 40}, and `rho - 1/3` in {1e-9, 1e-4, 0.01, 0.05, 0.15, 0.35, 0.667, 1.2,
2.0, 4.0, 8.0} — **704 points, no violations.** The largest value encountered
was `Q'' = -0.00184315`, at `n = 3`, `s = 1e-6`, `rho = 1/3 + 1e-9` — the
hardest corner of the grid on both counts, smallest `s` and `rho` closest to
`1/3`. It agrees with the small-`s` prediction to four digits:
`sqrt(1e-6) * L(3) = 1e-3 * (-1.84303) = -0.00184303`. So the grid maximum is
not an accident of sampling; it is the `s -> 0` limit above, and Target 1 is
what controls it.

Derivatives are **not** finite-differenced. `u', u'', v', v''` come from
`d/db E_b f = Cov_b(T,f)` and
`d^2/db^2 E_b f = E_b[((T - mu_b)^2 - sigma_b^2) f]`, so each is an exact
truncated-normal moment ratio evaluated by quadrature. This matters: at
`rho = 1/3 + 1e-9` a finite-difference `Q''` would be dominated by
differencing error.

**Status: the concavity target is unrefuted and looks right.** It is not
proved. The scan is evidence, not a proof, and no scan can exclude a violation
in an unsampled corner.

---

## 2. Conjecture 16: an asymptotic version is provable now

Let `n = m - 1` and `Qb` denote the standard normal survival function.

> **Theorem.** Fix `m >= 4` and `x > 0`. Let `M_m(x) = sup_G F_x(G)`, the
> supremum of `F_x(G) = sup_{E>0} A_x(EG)` over **every** Gram matrix with unit
> diagonal. Then
>
> ```
> 0  <=  M_m(x) - F_x(G_tri)  <=  binom(n,2) * Qb( sqrt(8x/3) ).
> ```

**Proof.** *Upper bound.* For any configuration and any `E > 0`, a union bound
in (1.3) gives `A_x(EG) <= (1/m) sum_i sum_{j != i} P{l_ij > x}`. If
`s_ij = 0` the competitor is identical to the true signal, `l_ij == 0`, and it
contributes nothing for `x > 0`. If `s_ij > 0` then `l_ij ~ N(-s_ij, 2 s_ij)`
and

```
P{l_ij > x} = Qb( (x + s_ij) / sqrt(2 s_ij) )  <=  Qb( sqrt(2x) ),
```

because `x + s >= 2 sqrt(x s)` by AM-GM, so
`(x + s)/sqrt(2s) >= 2 sqrt(xs)/sqrt(2s) = sqrt(2x)`. Each of the `m` inner
sums has at most `n` terms, so `A_x(EG) <= n Qb(sqrt(2x))` for every `E`, hence
`M_m(x) <= n Qb(sqrt(2x))`.

*Lower bound.* Take `E = E* = x(m-1)/m`, so `s = x`. At `G_tri` every
`s_ij = s = x`, and the individual tail attains the bound above with equality:
`Qb((x + x)/sqrt(2x)) = Qb(sqrt(2x))`. For `j != k`, both distinct from `i`,

```
Cov(l_ij, l_ik) = <mu_j - mu_i, mu_k - mu_i>
                = ( ||mu_j - mu_i||^2 + ||mu_k - mu_i||^2 - ||mu_j - mu_k||^2 ) / 2
                = (2s + 2s - 2s)/2 = s,
```

while `Var(l_ij) = 2s`, so the standardized scores `X_j = (l_ij + s)/sqrt(2s)`
are standard normal with correlation `1/2` and the exceedance threshold is
`t = (x + s)/sqrt(2s) = sqrt(2x)`. Then `Var(X_j + X_k) = 2 + 2(1/2) = 3` and

```
P{X_j > t, X_k > t}  <=  P{X_j + X_k > 2t}  =  Qb( 2t/sqrt(3) )  =  Qb( sqrt(8x/3) ).
```

Bonferroni gives
`F_x(G_tri) >= A_x(E* G_tri) >= n Qb(sqrt(2x)) - binom(n,2) Qb(sqrt(8x/3))`.
Subtracting from the upper bound proves the right inequality; the left holds
because `G_tri` is admissible. `[]`

Two remarks the proof sketch in the plan does not make explicit.

**The upper bound never used centering or cyclicity.** It only used
`s_ij >= 0`. So `M_m(x)` may be taken over all unit-diagonal Grams — the whole
slice `G^0_m` and more — not just the cyclic family. The conclusion is
correspondingly stronger than Conjecture 16's own scope.

**The bound is vacuous unless `x` is large compared to `log n`.** Since

```
Qb(sqrt(8x/3)) / Qb(sqrt(2x))  ->  (sqrt(3)/2) e^{-x/3},
```

the lower bound is positive only once `binom(n,2)(sqrt(3)/2)e^{-x/3} < n`,
i.e. roughly `x > 3 log n`. Below that threshold the theorem says nothing. For
fixed `m` this is harmless and gives

```
F_x(G_tri) / M_m(x)  =  1 - O(e^{-x/3}),
```

with explicit constant: the relative gap is at most
`((n-1) sqrt(3)/4) e^{-x/3} (1 + o(1))`. But it is exactly why this is a
fixed-`m` statement, and why it does not begin to address the `m -> infinity`
question the paper lists as open.

The bracket was checked numerically at `m = 4, 5` for
`x` in {0.5, 1, 2, 4, 8}: `F_x(G_tri)` lies inside it in every case, and the
profiling energy comes out at `E* = x(m-1)/m` to 6 digits, confirming Lemma 7
numerically. (At `x = 16` the bracket is narrower than the `1e-9` absolute
accuracy of the orthant integrator, so that row is not a test of anything.)

---

## 3. Exact global maximality at `m = 4`: numerics, and a trap

For `m = 4` the cyclic profile is `w = (w, 1-w)` with `G_{r=1} = -(1-w)` and
`G_{r=2} = -w + (1-w)`; the simplex is `w* = 2/3`. Sweeping `w` on a 0.05 grid
and profiling over `E`:

| x | F(w*) | best off-simplex | at w | E* found | x(m-1)/m |
|---|---|---|---|---|---|
| 0.5 | 0.322220605 | 0.322087719 | 0.65 | 0.375000 | 0.375 |
| 1.0 | 0.177207099 | 0.177176267 | 0.65 | 0.750000 | 0.75 |
| 2.0 | 0.057466558 | 0.057449731 | 0.65 | 1.500000 | 1.5 |
| 4.0 | 0.006529921 | 0.006526833 | 0.65 | 3.000000 | 3.0 |

`w* = 2/3` wins at every margin tested, and `F_x` looks unimodal in `w`. So
Conjecture 16 is numerically supported at `m = 4`, and the proposed route —
`F_x' > 0` below `2/3`, `F_x' < 0` above — is consistent with the data on the
interior.

**But the boundary is not a limit of the interior computation.** At `w = 1`
the Gram matrix has rank 2: the signals are `(1,0), (0,1), (-1,0), (0,-1)`, the
biorthogonal set in `R^2`, and the competitor covariance
`E * [[2,2,0],[2,4,2],[0,2,2]]` is singular. A Cholesky-based orthant routine
with a small regularizer returns nonsense there, and in a first pass it
produced an apparent *upward jump* at `w = 1` that looked exactly like a second
local maximum at the boundary — at `x = 4`, a spurious `0.00607` against
`0.00555` at `w = 0.99`. Computing that point from its own geometry instead
(with `W, V` iid `N(0,E)`, the event is
`W + |V| > x + E` or `W > x/2 + E`, a one-dimensional quadrature) gives the
true values:

| x | F(w=1), exact rank-2 | F(w*) |
|---|---|---|
| 2 | 0.0506007 | 0.0574666 |
| 4 | 0.0054964 | 0.0065299 |
| 6 | 0.0006284 | 0.0007723 |

No second local maximum; the simplex still wins. The lesson for the proof plan
is concrete: a `w`-derivative argument on the cyclic family has to treat the
rank-deficient boundary profiles separately, because the score covariance
degenerates there and the profiled functional is not accessible by the same
representation used on the interior. The plan's "appropriate treatment of
optimizer switches and boundary limits" is not a formality — it is where the
numerics break first.

---

## 4. Where this leaves the sequence

1. **Target 1** (`(n-1)u(0)/2 + n(n+1)v'(0) < 0`) is the smallest true
   sub-problem of the concavity route, it isolates the one unconditionally
   positive term in `Q''`, and identity (B) plus `R' = -R(T+R)` is the tool for
   it. Doing `n = 3` first is cheap.
2. The asymptotic global bound of section 2 is proved and can go into the paper
   as stated, with the two strengthenings noted (all unit-diagonal Grams; the
   explicit `x > 3 log n` threshold).
3. Exact `m = 4` global maximality: numerics support it; the boundary needs its
   own argument.

Reproduce with:

```bash
python scripts/curvature_conjectures.py --quick     # a few minutes
python scripts/curvature_conjectures.py             # full 704-point scan, ~20 min
```
