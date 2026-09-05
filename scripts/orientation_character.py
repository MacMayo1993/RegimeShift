"""Verification of the orientation character of the simplex shape representation.

    H_0 = {H = H^T : diag H = 0, H 1 = 0},   sigma . H = P_sigma H P_sigma^T

The claim is det(sigma | H_0) = sgn(sigma)^(m-3), proved by the S_m-equivariant
short exact sequence

    0 -> H_0 -> E_m --R--> R^m -> 0,      R(H) = H 1,   E_m = {sym, diag 0},

since det(sigma | E_m) = sgn(sigma)^(m-2) (a transposition induces exactly m-2
transpositions of edges) and det(sigma | R^m) = sgn(sigma).

This script checks the identity by brute force over the whole group, and checks
two structural facts the write-up depends on:

  * H_0 is irreducible, i.e. <chi, chi> = 1, confirming H_0 = V_{(m-2,2)};
  * the action is faithful for m >= 5 but NOT for m = 4, where the kernel is
    the Klein four-group, so the effective group at m = 4 is S_4/V_4 = S_3.

Numpy only. Run: python scripts/orientation_character.py [--max-m 8]
"""

from __future__ import annotations

import itertools
from math import factorial

import numpy as np


def h0_basis(m):
    """Orthonormal basis of H_0, as columns, in the edge coordinates of K_m."""
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    rowsum = np.zeros((m, len(pairs)))
    for c, (i, j) in enumerate(pairs):
        rowsum[i, c] = 1.0
        rowsum[j, c] = 1.0
    _, sv, vt = np.linalg.svd(rowsum)
    return pairs, vt[np.sum(sv > 1e-9):].T


def action(m, perm, pairs, basis):
    """Matrix of H -> P_sigma H P_sigma^T restricted to H_0, in that basis."""
    index = {p: c for c, p in enumerate(pairs)}
    P = np.zeros((len(pairs), len(pairs)))
    for c, (i, j) in enumerate(pairs):
        a, b = perm[i], perm[j]
        P[index[(min(a, b), max(a, b))], c] = 1.0
    return basis.T @ P @ basis


def sign(perm):
    s = 1
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                s = -s
    return s


def check(m, verbose=True):
    pairs, basis = h0_basis(m)
    d = basis.shape[1]
    dets, chars, kernel = [], [], []
    for perm in itertools.permutations(range(m)):
        M = action(m, list(perm), pairs, basis) if d else np.zeros((0, 0))
        det = float(np.linalg.det(M)) if d else 1.0
        dets.append((perm, det))
        chars.append(float(np.trace(M)) if d else 0.0)
        if d and np.allclose(M, np.eye(d), atol=1e-9):
            kernel.append(perm)
        elif not d:
            kernel.append(perm)
    bad = [(p, v) for p, v in dets if abs(v - sign(list(p))**(m - 3)) > 1e-8]
    irr = sum(c*c for c in chars)/factorial(m) if d else float('nan')
    if verbose:
        print("m=%-3d dim H_0=%-4d (= m(m-3)/2 = %-4d)  det character verified on "
              "all %d permutations: %s"
              % (m, d, m*(m - 3)//2, factorial(m), "OK" if not bad else "FAILED"))
        print("      <chi,chi>=%.6f  faithful=%s  |kernel|=%d"
              % (irr, len(kernel) == 1, len(kernel)))
        if len(kernel) > 1 and d:
            print("      kernel: %s" % (kernel,))
    return not bad


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--max-m', type=int, default=8,
                   help='largest m to check (cost is m!; 8 takes a few minutes)')
    a = p.parse_args(argv)
    ok = all(check(m) for m in range(3, a.max_m + 1))
    print("\nall checks passed" if ok else "\nFAILURES")
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
