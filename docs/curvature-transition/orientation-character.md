# The orientation character of the simplex shape representation

Companion to `conjecture-notes.md`. This one is not about the conjectures; it
records a symmetry fact about the tangent space `H_0` and, more importantly,
separates it cleanly from the margin-driven curvature transition of Theorem 3.

Verified by `scripts/orientation_character.py`.

---

## 1. The proposition, and its verification

> **Proposition.** Let `H_0 = {H = H^T : diag H = 0, H 1 = 0}` carry the
> relabeling action `sigma . H = P_sigma H P_sigma^T`. For `m >= 3`,
>
> ```
> det( sigma | H_0 ) = sgn(sigma)^(m-3).
> ```

The proof by short exact sequence is correct and needs nothing added:

```
0 -> H_0 -> E_m --R--> R^m -> 0,    E_m = {sym, diag 0},   R(H) = H1,
```

is `S_m`-equivariant because `P_sigma^T 1 = 1`; `R` is onto for `m >= 3` since
`R(e_ij) = e_i + e_j` and `e_i = [(e_i+e_j) + (e_i+e_k) - (e_j+e_k)]/2`; `E_m`
is the permutation representation on edges, where a transposition induces
exactly the `m-2` swaps `{1,k} <-> {2,k}`, `k = 3..m`, giving
`det(sigma|E_m) = sgn(sigma)^(m-2)`; `det(sigma|R^m) = sgn(sigma)`;
multiplicativity along the sequence and `sgn^{-1} = sgn` finish it. The
`m >= 3` hypothesis is needed and worth writing down — surjectivity of `R`
uses three distinct indices, and at `m = 2` it genuinely fails.

Checked by brute force on the **entire group** for `m = 3..8`
(`det(sigma|H_0)` against `sgn(sigma)^(m-3)` for all `sigma`, 40320
permutations at `m = 8`): no discrepancies, and `dim H_0 = m(m-3)/2` in every
case. Two further facts came out of the same computation:

| m | dim H_0 | `<chi,chi>` | action faithful? | kernel |
|---|---|---|---|---|
| 4 | 2 | 1 | **no** | `V_4 = {e, (12)(34), (13)(24), (14)(23)}` |
| 5 | 5 | 1 | yes | trivial |
| 6 | 9 | 1 | yes | trivial |
| 7 | 14 | 1 | yes | trivial |

`<chi,chi> = 1` confirms `H_0` is irreducible, i.e. `H_0 = V_{(m-2,2)}`
independently of any character table — so the determinant computation is
indeed an orientation character for that Specht module.

---

## 2. Three places the statement needs tightening

Everything below is a correction of scope, not of the algebra.

### (a) At `m = 4` the action is not faithful

The kernel is the Klein four-group, so the effective group acting on `H_0` is
`S_4/V_4 = S_3`, realized on `R^2` as its 2-dimensional reflection
representation. The determinant identity is unaffected — `(12)(34)` is even and
in the kernel, and `det = +1` both ways — but the *orbifold* sentence has to be
read through the image: the linearized quotient at the simplex is `R^2/S_3`, a
reflection orbifold, not `R^2/S_4`. The orientation-reversing conclusion
survives, since a transposition still maps to a reflection.

### (b) The step from a character to *the* orientation bundle needs the slice to be affine

Section 8 of the argument forms `L_or = U x_{S_m} det(H_0)` and calls it the
orientation local system of `U/S_m`. As written that is a flat line bundle
associated to a character — which is not by itself the orientation bundle,
because `H_0` is the tangent space at `G_tri` specifically, and `G_tri` is not
in `U` (the action is not free there).

What closes the gap is that the slice is **affine**:
`{G : diag G = 1, G1 = 0}` is an affine subspace with direction `H_0`, and
`G^0_m` is its intersection with the PSD cone. So at every relative-interior
point the tangent space is canonically `H_0`, and the derivative of the
`sigma`-action at *any* such point is exactly `rho(sigma)` — the same
representation, not merely a conjugate of it. Hence for open `S_m`-invariant
`U` in the relative interior on which the action is free, the orientation
bundle of `U/S_m` really is `U x_{S_m} det(H_0)`, with monodromy
`sgn^{m-3}`. With that sentence added, section 8 is correct.

The relative-interior restriction is not a technicality. On the rank-deficient
boundary `G^0_m` is not a manifold, and there is no tangent space to orient —
the same boundary degeneracy that, at `m = 4`, produced the spurious second
local maximum documented in `conjecture-notes.md` section 3.

### (c) Section 11 proves less than the available argument

The argument given restricts to an interval `I` in `rho` free of zeros of
`I_m`, notes `pi_1(I) = 0`, and concludes `w_1(E) = 0`. That is true but
vacuous: *every* bundle over an interval is trivial, so this would "prove" the
same thing for a Hessian with genuinely twisting eigenspaces.

The real content is stronger and does not need a simply connected parameter
domain at all. By Theorem 1 the Hessian, transported through the Frobenius
metric, is

```
L_{s,x} = c(s,x) * Id_{H_0},        c(s,x) = -I_m(s,x)/(8 m s^2),
```

a **scalar** operator, because `||H||_F^2` carries no dependence on `s` or `x`.
Every vector of `H_0` is an eigenvector at every parameter value. So the
eigenbundle is the product `H_0 x P` over *any* parameter space `P`, however
complicated its topology, and there is no eigenvector monodromy to speak of.
The curvature reversal is a crossing of the discriminant `{I_m = 0}`, where the
whole form vanishes simultaneously — not transport around it. That is the
statement to make; the `pi_1(I) = 0` version should be dropped.

Relatedly, the parity table offered as a "falsification test" is corroboration,
not proof. That the curvature transition exists for both parities while the
orientation character depends on parity is good evidence the two are
unrelated, but the actual argument is the scalar form above, which settles it
outright. Worth keeping the table, worth not calling it a mechanism-level
proof.

---

## 3. What is safe to state

Both of these are fully supported:

> **Proposition.** Under `H -> P_sigma H P_sigma^T`,
> `det(sigma|H_0) = sgn(sigma)^(m-3)` for `m >= 3`. Hence for odd `m` every
> relabeling preserves orientation of `H_0`; for even `m` odd relabelings
> reverse it. On an open `S_m`-invariant subset of the relative interior of
> `G^0_m` where the action is free, the orientation bundle of the quotient is
> the flat line bundle with monodromy `sgn^{m-3}`. At `m = 4` the action
> factors through `S_4/V_4 = S_3`.

> **Remark.** This parity phenomenon is independent of the curvature
> transition. The intrinsic shape Hessian is `c(s,x) ||.||_F^2` on `H_0`, a
> scalar multiple of the identity at every parameter value, so its eigenspaces
> are constant and cannot twist with the margin; only the scalar
> `c = -I_m/(8ms^2)` changes sign, and it does so by passing through zero.

What should stay conjectural is the global claim: `det rho != 1` does not give
`w_1(U/S_m) != 0` for any particular `U`. That needs an actual loop in
`U/S_m` whose monodromy is an odd permutation, i.e. an odd element in the image
of `pi_1(U/S_m) -> S_m`, and no such loop is exhibited. The local orbifold
statement at the simplex needs no such input, because an orientation-reversing
element sits in the isotropy group directly.

Two smaller notes. `m(m-3)/2 = 3` has no integer root (`m = (3 + sqrt(33))/2`),
so `P(H_0) = RP^{m(m-3)/2 - 1}` is never `RP^2` — the observation is right, and
`RP^1, RP^4, RP^8` for `m = 4, 5, 6`. And for the representation theory of
`S_m` the natural citations are the standard references (Fulton-Harris, or
James's *The Representation Theory of the Symmetric Groups*) rather than a
search-aggregator link; the proof above avoids needing them at all, which is
its main virtue.
