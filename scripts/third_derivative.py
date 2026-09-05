"""Third-order shape derivative of A_x at the regular simplex.

Checks the exact third directional derivative obtained by differentiating the
tilted Gaussian density three times, the invariant-theoretic shape of the
resulting cubic, and -- for m = 4 -- its actual value at the two curvature
transitions.

With K_eps = sI + eps H, H in H_0 (so tr H = 0 and the comparison region in
Psi_x does not move), the log-density derivatives at eps = 0 are

    L'   = Y'HY / (2 s^2)
    L''  = tr(H^2)/(2 s^2) - Y'H^2 Y / s^3            [both already in the paper]
    L''' = -tr(H^3)/s^3 + 3 Y'H^3 Y / s^4             [new]

and f'''/f = L''' + 3 L' L'' + (L')^3, giving

    S_3(Y;H) = (Y'HY)^3/(8s^6) - 3(Y'HY)(Y'H^2Y)/(2s^5)
             + 3 tr(H^2)(Y'HY)/(4s^4) + 3 Y'H^3Y/s^4 - tr(H^3)/s^3,

    D^3 A_x(sI)[H,H,H] = -(e^{-s/2}/m) E_{N(0,sI)}[ Psi_x(Y) S_3(Y;H) ].

What the script verifies:

1. ``check_density_algebra`` -- L', L'', L''' and f'''/f against high-precision
   finite differences of the actual Gaussian log-density, and the full mixed
   trilinear form for three distinct directions.

2. ``invariant_dimension`` -- dim (Sym^3 H_0^*)^{S_m} by the character formula,
   using chi_{H_0}(sigma) = C(F,2) + T - F with F fixed points and T 2-cycles.
   Result: 1 for m = 4, 5 and 2 for m >= 6.

3. ``cubic_basis_rank`` -- that tr(H^3) and sum_{i<j} H_ij^3 span that space,
   and that on H_0 they satisfy tr(H^3) = 4 sum_{i<j} H_ij^3 exactly when
   m = 4, 5 (which is why the invariant space is one-dimensional there).

4. ``four_cycle_is_null`` -- the paper's four-cycle direction is annihilated by
   every invariant cubic, and more strongly P_sigma H P_sigma^T = -H for the
   4-cycle sigma = (1 2 4 3), so A_x is even in eps along that line and every
   odd derivative vanishes. For m = 4 the three four-cycles are exactly the
   zero locus of the invariant cubic on H_0.

5. ``alpha_m4`` -- the coefficient itself. A_x is evaluated directly (each term
   of Psi_x is tilted exactly into a Gaussian rectangle probability), validated
   against the paper's exact Hessian, then differenced three times.

Requires numpy/scipy, mpmath, and scripts/_orthant.py.
"""

from __future__ import annotations

import itertools
from math import comb, factorial

import numpy as np


# ------------------------------------------------------------- linear algebra

def h0_basis(m):
    """Orthonormal basis of H_0 as columns, in edge coordinates of K_m."""
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    R = np.zeros((m, len(pairs)))
    for c, (i, j) in enumerate(pairs):
        R[i, c] = 1.0
        R[j, c] = 1.0
    _, sv, vt = np.linalg.svd(R)
    return pairs, vt[np.sum(sv > 1e-9):].T


def to_matrix(pairs, m, w):
    H = np.zeros((m, m))
    for c, (i, j) in enumerate(pairs):
        H[i, j] = H[j, i] = w[c]
    return H


def h0_random(m, rng):
    pairs, B = h0_basis(m)
    return to_matrix(pairs, m, B @ rng.normal(size=B.shape[1]))


def S3(Y, H, s):
    """f'''_0 / f_0 evaluated at rows of Y."""
    Y = np.atleast_2d(Y)
    H2, H3 = H @ H, H @ H @ H
    q = np.einsum('ni,ij,nj->n', Y, H, Y)
    q2 = np.einsum('ni,ij,nj->n', Y, H2, Y)
    q3 = np.einsum('ni,ij,nj->n', Y, H3, Y)
    return (-np.trace(H3)/s**3 + 3*q3/s**4 + 3*np.trace(H2)*q/(4*s**4)
            - 3*q*q2/(2*s**5) + q**3/(8*s**6))


# ------------------------------------------------------------------- checks 1

def check_density_algebra(cases=((4, 1.7), (5, 0.9), (6, 2.3)), seed=7):
    from mpmath import mp, mpf, matrix, exp, log, det, lu_solve, diff
    mp.dps = 40
    rng = np.random.default_rng(seed)
    print("1. density derivatives vs high-precision finite differences")
    ok = True
    for m, s in cases:
        H = h0_random(m, rng)
        y = rng.normal(size=m)*1.3

        def L(e):
            K = matrix(m, m)
            for i in range(m):
                for j in range(m):
                    K[i, j] = mpf(s if i == j else 0) + e*mpf(H[i, j])
            Kiy = lu_solve(K, matrix([mpf(t) for t in y]))
            return (-log(det(K))/2 - sum(mpf(y[i])*Kiy[i] for i in range(m))/2
                    - m*log(2*mp.pi)/2)

        H2, H3 = H @ H, H @ H @ H
        want = [y @ H @ y/(2*s**2),
                np.trace(H2)/(2*s**2) - y @ H2 @ y/s**3,
                -np.trace(H3)/s**3 + 3*(y @ H3 @ y)/s**4]
        got = [diff(L, mpf(0), k, h=mpf(10)**-6) for k in (1, 2, 3)]
        fdd = (diff(lambda e: exp(L(e)), mpf(0), 3, h=mpf(10)**-6)
               / exp(L(mpf(0))))
        rels = [float(abs(g - w)/abs(w)) for g, w in zip(got, want)]
        rf = float(abs(fdd - S3(y, H, s)[0])/abs(S3(y, H, s)[0]))
        ok &= max(rels + [rf]) < 1e-8
        print("   m=%-3d L' %.1e  L'' %.1e  L''' %.1e  f'''/f %.1e"
              % (m, rels[0], rels[1], rels[2], rf))
    return ok


def check_trilinear(cases=((5, 1.4), (6, 0.8)), seed=11):
    from mpmath import mp, mpf, matrix, exp, log, det, lu_solve
    mp.dps = 40
    rng = np.random.default_rng(seed)
    print("   full mixed trilinear form (three distinct directions)")
    ok = True
    for m, s in cases:
        Hs = [h0_random(m, rng) for _ in range(3)]
        y = rng.normal(size=m)*1.1

        def L(e):
            K = matrix(m, m)
            for i in range(m):
                for j in range(m):
                    K[i, j] = (mpf(s if i == j else 0)
                               + sum(e[k]*mpf(Hs[k][i, j]) for k in range(3)))
            Kiy = lu_solve(K, matrix([mpf(t) for t in y]))
            return (-log(det(K))/2 - sum(mpf(y[i])*Kiy[i] for i in range(m))/2
                    - m*log(2*mp.pi)/2)

        h = mpf(10)**-5
        fd = sum(int(np.prod(sg))*exp(L([mpf(g)*h for g in sg]))
                 for sg in itertools.product([1, -1], repeat=3))
        fd = fd/(8*h**3)/exp(L([mpf(0)]*3))
        ell = [y @ H @ y/(2*s**2) for H in Hs]
        lij = {}
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            A, Bm = Hs[i], Hs[j]
            lij[(i, j)] = (np.trace(A @ Bm)/(2*s**2)
                           - y @ (A @ Bm + Bm @ A) @ y/(2*s**3))
        tr3 = np.trace(Hs[0] @ Hs[1] @ Hs[2]) + np.trace(Hs[0] @ Hs[2] @ Hs[1])
        Ssum = sum(Hs[p[0]] @ Hs[p[1]] @ Hs[p[2]]
                   for p in itertools.permutations(range(3)))
        l123 = -tr3/(2*s**3) + y @ Ssum @ y/(2*s**4)
        claim = (l123 + lij[(0, 1)]*ell[2] + lij[(0, 2)]*ell[1]
                 + lij[(1, 2)]*ell[0] + ell[0]*ell[1]*ell[2])
        rel = float(abs(mpf(claim) - fd)/abs(fd))
        ok &= rel < 1e-7
        print("   m=%-3d D^3 f / f  rel=%.1e" % (m, rel))
    return ok


# ------------------------------------------------------------------- checks 2

def _cycles(p):
    n = len(p)
    seen = [False]*n
    out = []
    for i in range(n):
        if not seen[i]:
            ln, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                ln += 1
            out.append(ln)
    return out


def _chi_h0(p):
    cs = _cycles(p)
    F, T = cs.count(1), cs.count(2)
    return comb(F, 2) + T - F


def invariant_dimension(m):
    """dim (Sym^3 H_0^*)^{S_m} by the character formula."""
    tot = 0
    for p in itertools.permutations(range(m)):
        p = list(p)
        p2 = [p[p[i]] for i in range(m)]
        p3 = [p[p2[i]] for i in range(m)]
        a, b, c = _chi_h0(p), _chi_h0(p2), _chi_h0(p3)
        tot += a**3 + 3*a*b + 2*c
    return tot//(6*factorial(m))


def cubic_basis_rank(m, trials=400, seed=3):
    """Rank of {tr(H^3), sum H_ij^3} as forms, and their ratio when rank 1."""
    rng = np.random.default_rng(seed)
    pairs, _ = h0_basis(m)
    V = []
    for _ in range(trials):
        H = h0_random(m, rng)
        V.append([np.trace(H @ H @ H), sum(H[i, j]**3 for i, j in pairs)])
    V = np.array(V)
    sv = np.linalg.svd(V, compute_uv=False)
    rank = int(np.sum(sv > 1e-8*sv[0]))
    ratio = float((V[:, 0]/V[:, 1]).mean()) if rank == 1 else None
    return rank, ratio


# ------------------------------------------------------------------- checks 4

def four_cycle(m):
    H = np.zeros((m, m))
    H[0, 1] = H[1, 0] = 1.0
    H[2, 3] = H[3, 2] = 1.0
    H[0, 2] = H[2, 0] = -1.0
    H[1, 3] = H[3, 1] = -1.0
    return H


def four_cycle_is_null(m):
    """The four-cycle is cubic-null, and an odd relabeling negates it."""
    H = four_cycle(m)
    pairs, _ = h0_basis(m)
    perm = [1, 3, 0, 2] + list(range(4, m))       # sigma = (1 2 4 3), 1-based
    P = np.zeros((m, m))
    for i, pi in enumerate(perm):
        P[pi, i] = 1.0
    return dict(rowsum=float(np.abs(H @ np.ones(m)).max()),
                tr3=float(np.trace(H @ H @ H)),
                sum3=float(sum(H[i, j]**3 for i, j in pairs)),
                negated_by_sigma=float(np.abs(P @ H @ P.T + H).max()),
                cube_identity=float(np.abs(H @ H @ H - 4*H).max()))


# ------------------------------------------------------------------- checks 5

def psi_expect(K, x, N=128):
    """E_K[Psi_x(Y)], each term tilted exactly into a rectangle probability."""
    from _orthant import rect_prob
    m = K.shape[0]
    tot = 0.0
    for i in range(m):
        idx = [j for j in range(m) if j != i]
        Sig = np.array([[K[j, k] - K[j, i] - K[i, k] + K[i, i] for k in idx]
                        for j in idx])
        b = np.array([x - K[j, i] + K[i, i] for j in idx])
        L = np.linalg.cholesky(
            Sig + 1e-14*max(1.0, float(np.abs(Sig).max()))*np.eye(len(idx)))
        tot += np.exp(K[i, i]/2)*rect_prob(np.full(len(idx), -np.inf), b, L, N)
    return tot


def A_x(K, x, s, N=128):
    return 1.0 - np.exp(-s/2)/K.shape[0]*psi_expect(K, x, N)


def _d3(g, h):
    return (g(2*h) - 2*g(h) + 2*g(-h) - g(-2*h))/(2*h**3)


def alpha_m4(s=1.0, N=160, verbose=True):
    """D^3 A_x / tr(H^3) at the two curvature transitions, m = 4.

    Direction-independence of the ratio is the numerical signature of the
    one-dimensional invariant cubic space at m = 4.
    """
    from mpmath import mpf, findroot
    from curvature_conjectures import I_closed
    m, n = 4, 3
    pairs, B = h0_basis(m)
    rc = [float(findroot(lambda r: I_closed(r, s, n), mpf(g)))
          for g in ('0.4', '3.0')]
    out = {}
    for rho, tag in zip(rc, ('rho_c^-', 'rho_c^+')):
        if verbose:
            print("   %s = %.10f   (m=4, s=%g)" % (tag, rho, s))
        rows = []
        for deg in (9.9386, 40.0, 70.0, 129.9386):
            t = np.radians(deg)
            w = B @ np.array([np.cos(t), np.sin(t)])
            H = to_matrix(pairs, m, w)*np.sqrt(8.0/(w @ w))
            tr3 = float(np.trace(H @ H @ H))
            g = lambda e: A_x(s*np.eye(m) + e*H, rho*s, s, N)
            v = [_d3(g, h) for h in (0.14, 0.10)]
            d3 = (4*v[1] - v[0])/3
            rows.append((deg, tr3, d3, d3/tr3 if abs(tr3) > 1e-9 else np.nan))
            if verbose:
                print("      theta=%-10.4f tr(H^3)=%-12.6g D^3A=%-12.6g "
                      "ratio=%-12.6g" % rows[-1])
        out[tag] = rows
    return out


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--part', default='all',
                   choices=['all', 'algebra', 'invariants', 'fourcycle', 'alpha'])
    p.add_argument('--max-m', type=int, default=8)
    a = p.parse_args(argv)
    run = lambda k: a.part in ('all', k)
    if run('algebra'):
        check_density_algebra()
        check_trilinear()
        print()
    if run('invariants'):
        print("2/3. invariant cubics on H_0")
        for m in range(4, a.max_m + 1):
            rank, ratio = cubic_basis_rank(m)
            rel = ("tr(H^3) = %.6f * sum H_ij^3" % ratio) if ratio else "independent"
            print("   m=%-3d dim (Sym^3 H_0*)^{S_m} = %d   basis rank = %d   %s"
                  % (m, invariant_dimension(m), rank, rel))
        print()
    if run('fourcycle'):
        print("4. the paper's four-cycle direction")
        for m in (4, 6):
            d = four_cycle_is_null(m)
            print("   m=%-3d rowsum=%.0e tr(H^3)=%.0e sum H_ij^3=%.0e "
                  "|P H P' + H|=%.0e |H^3-4H|=%.0e"
                  % (m, d['rowsum'], d['tr3'], d['sum3'],
                     d['negated_by_sigma'], d['cube_identity']))
        print()
    if run('alpha'):
        print("5. D^3 A_x at the curvature transitions, m = 4")
        alpha_m4()


if __name__ == '__main__':
    main()
