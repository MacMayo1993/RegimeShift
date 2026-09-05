"""Fourth variation of A_x along the four-cycle direction.

Along the paper's direction (H_12 = H_34 = 1, H_13 = H_24 = -1) the whole
covariance perturbation collapses to a single correlation. With
u = e_1 - e_4 and v = e_2 - e_3 one has H = u v' + v u', u'v = 0,
|u|^2 = |v|^2 = 2, so for K_eps = sI + eps H the standardized pair

    U = (Y_1 - Y_4)/sqrt(2s),   V = (Y_2 - Y_3)/sqrt(2s)

is bivariate normal with correlation r = 2 eps / s, and everything orthogonal
to span{u,v} is untouched. Mehler's expansion then gives

    d^k/d eps^k (f_eps / f_0) |_0 = (2/s)^k He_k(U) He_k(V),

so with P = (Y_1-Y_4)^2 and Q = (Y_2-Y_3)^2,

    k=2:  (P - 2s)(Q - 2s) / s^4                       [the paper's Hessian]
    k=4:  (P^2-12sP+12s^2)(Q^2-12sQ+12s^2) / s^8       [this module]

Conditioning as in Section 5.4 (tilt the base coordinate, standardize, truncate
the others at t) reduces the k=4 case to a one-dimensional integral

    K_m(s,x) = int W_b(t) [ (m-4) D(t)^2 + 4 C(t) D(t) ] dt,
    W_b(t)   = phi(t-b) Phi(t)^n,   b = sqrt(s) + x/sqrt(s),   n = m-1,

the exact quartic analogue of the paper's I_m = int W_b[(m-4)B^2 + 4AB], with

    D(t) = -2 lam [ t^3 - 3t + lam (t^2 - 4) ]
    C(t) = He_4(c) + lam [ 4c^3 - 6c^2 t + 4c t^2 - 16c - t^3 + 9t ],
    lam = phi/Phi,  c = t - x/sqrt(s),  He_4(c) = c^4 - 6c^2 + 3.

Then  D^4 A_x(sI)[H,H,H,H] = -K_m(s,x) / (m s^4).

D(t) is 4 E[He_4((X-X')/sqrt 2)] for X, X' iid N(0,1) truncated at t (both
vertices of a pair ordinary); C(t) is the same with one vertex the transmitted
base. Four of the m bases lie in {1,2,3,4} and give C*D, the other m-4 give
D^2 -- exactly the counting behind the Hessian formula.

Because a transposition sigma = (2 3) satisfies P_sigma H P_sigma' = -H and
A_x is S_m-invariant, eps -> A_x(sI + eps H) is exactly even, so every odd
derivative vanishes along this line and the quartic is the first correction
once the Hessian vanishes -- on this slice.

Run: python scripts/quartic_four_cycle.py --help
"""

from __future__ import annotations

import numpy as np
from mpmath import mp, mpf, exp, log, sqrt, quad, erfc

mp.dps = 30


def _logphi(t): return -t*t/2 - log(2*mp.pi)/2
def _Phi(t):    return erfc(-t/sqrt(2))/2
def _lam(t):    return exp(_logphi(t))/_Phi(t)


def D_factor(t):
    """Ordinary-ordinary pair factor."""
    l = _lam(t)
    return -2*l*(t**3 - 3*t + l*(t*t - 4))


def C_factor(t, h):
    """Base-ordinary pair factor; h = x/sqrt(s), c = t - h."""
    l = _lam(t)
    c = t - h
    return ((c**4 - 6*c*c + 3)
            + l*(4*c**3 - 6*c*c*t + 4*c*t*t - 16*c - t**3 + 9*t))


def K_m(s, x, m, width=40):
    s, x = mpf(s), mpf(x)
    n = m - 1
    b = sqrt(s) + x/sqrt(s)
    h = x/sqrt(s)
    f = lambda t: (exp(-(t-b)**2/2 - log(2*mp.pi)/2 + n*log(_Phi(t)))
                   * ((m-4)*D_factor(t)**2 + 4*C_factor(t, h)*D_factor(t)))
    return quad(f, [b-width, b-6, b-1, b+1, b+6, b+width])


def D4A(s, x, m):
    """D^4 A_x(sI)[H,H,H,H] along the four-cycle."""
    return -K_m(s, x, m)/(m*mpf(s)**4)


def four_cycle(m):
    H = np.zeros((m, m))
    H[0, 1] = H[1, 0] = 1.0
    H[2, 3] = H[3, 2] = 1.0
    H[0, 2] = H[2, 0] = -1.0
    H[1, 3] = H[3, 1] = -1.0
    return H


# ------------------------------------------------------------------- checks

def check_hermite(m=6, s=1.3, seed=5):
    """The k=2 and k=4 scores against finite differences of the true density."""
    from mpmath import matrix, det, lu_solve, diff
    mp.dps = 40
    H = four_cycle(m)
    y = np.random.default_rng(seed).normal(size=m)*1.2

    def L(e):
        K = matrix(m, m)
        for i in range(m):
            for j in range(m):
                K[i, j] = mpf(s if i == j else 0) + e*mpf(H[i, j])
        Kiy = lu_solve(K, matrix([mpf(t) for t in y]))
        return (-log(det(K))/2 - sum(mpf(y[i])*Kiy[i] for i in range(m))/2
                - m*log(2*mp.pi)/2)

    P, Q = (y[0]-y[3])**2, (y[1]-y[2])**2
    H2 = H @ H
    H4 = H2 @ H2
    out = {}
    out["tr"] = (float(np.trace(H2)), float(np.trace(H @ H2)), float(np.trace(H4)),
                 float(np.abs(H @ H2 - 4*H).max()), float(np.abs(H4 - 4*H2).max()))
    claim4 = 3*np.trace(H4)/s**4 - 12*(y @ H4 @ y)/s**5
    out["L4"] = float(abs(diff(L, mpf(0), 4, h=mpf(10)**-5) - mpf(claim4))/abs(claim4))
    for k, herm in ((2, (P-2*s)*(Q-2*s)/s**4),
                    (4, (P**2-12*s*P+12*s**2)*(Q**2-12*s*Q+12*s**2)/s**8)):
        fd = diff(lambda e: exp(L(e)), mpf(0), k, h=mpf(10)**-5)/exp(L(mpf(0)))
        out["f%d" % k] = float(abs(fd - mpf(herm))/abs(fd))
    mp.dps = 30
    return out


def check_against_finite_differences(cases=((4, 1.0, 0.6), (4, 1.0, 1.5),
                                            (5, 1.0, 1.0), (4, 2.0, 1.2),
                                            (6, 1.0, 0.8)), N=200):
    """-K_m/(m s^4) against direct 4th differences of A_x.

    A_x is evaluated by tilting each term of Psi_x into a Gaussian rectangle
    probability (see third_derivative.A_x). The 4th difference is extrapolated
    to h -> 0 with an h^2 + h^4 model on the three smallest steps; using large
    steps as well corrupts the fit, since the h^6 term is not negligible there.
    """
    from third_derivative import A_x
    rows = []
    for m, s, rho in cases:
        x = rho*s
        H = four_cycle(m)
        g = lambda e: A_x(s*np.eye(m) + e*H, x, s, N)
        d4 = lambda h: (g(2*h) - 4*g(h) + 6*g(0.0) - 4*g(-h) + g(-2*h))/h**4
        hs = [0.13*s, 0.10*s, 0.08*s]
        raw = [d4(h) for h in hs]
        a = np.linalg.solve(np.array([[1.0, h**2, h**4] for h in hs]),
                            np.array(raw))[0]
        th = float(D4A(s, x, m))
        rows.append((m, s, rho, th, float(a), abs(a-th)/abs(th)))
    return rows


def _bisect(f, lo, hi, iters=60):
    lo, hi = mpf(lo), mpf(hi)
    flo = f(lo)
    for _ in range(iters):
        mid = (lo + hi)/2
        fm = f(mid)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi)/2


def sign_law(ss=(0.25, 1.0, 4.0, 9.0), ms=(4, 5, 6, 10, 25)):
    """K_m > 0 at rho_c^- and K_m < 0 at rho_c^+ ?

    The second half is FALSE. See ``K_pieces`` for why: K_m has a manifestly
    positive part that grows linearly in m.
    """
    from curvature_conjectures import I_closed
    rows, ok = [], True
    for s in ss:
        for m in ms:
            n = m - 1
            F = lambda r: I_closed(r, s, n)
            r1 = _bisect(F, 0.34, 0.99)
            r2 = _bisect(F, 1.01, 12.0)
            k1, k2 = float(K_m(s, float(r1)*s, m)), float(K_m(s, float(r2)*s, m))
            good = k1 > 0 and k2 < 0
            ok &= good
            rows.append((m, s, float(r1), k1, float(r2), k2, good))
    return rows, ok


def K_pieces(s, x, m, width=40):
    """Split K_m = (m-4) * Idd + 4 * Icd,  Idd = int W_b D^2,  Icd = int W_b C D.

    Idd > 0 always, since D^2 >= 0 and D is not identically zero. So the first
    part is strictly positive and grows linearly in m, which is why
    K_m(rho_c^+) < 0 cannot survive large m.
    """
    s, x = mpf(s), mpf(x)
    n = m - 1
    b = sqrt(s) + x/sqrt(s)
    h = x/sqrt(s)
    W = lambda t: exp(-(t-b)**2/2 - log(2*mp.pi)/2 + n*log(_Phi(t)))
    pts = [b-width, b-6, b-1, b+1, b+6, b+width]
    return (quad(lambda t: W(t)*D_factor(t)**2, pts),
            quad(lambda t: W(t)*C_factor(t, h)*D_factor(t), pts))


def sign_law_breakdown(ss=(1.0,), ms=(6, 10, 12, 15, 25, 60)):
    """The two pieces of K_m at rho_c^+, showing where and why the law fails."""
    from curvature_conjectures import I_closed
    rows = []
    for s in ss:
        for m in ms:
            r2 = _bisect(lambda r: I_closed(r, s, m-1), 1.01, 14.0)
            idd, icd = K_pieces(s, float(r2)*s, m)
            rows.append((m, s, float(r2), float(idd), float(icd),
                         float((m-4)*idd + 4*icd)))
    return rows


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--part', default='all',
                   choices=['all', 'hermite', 'finite-diff', 'table', 'signs'])
    a = p.parse_args(argv)
    run = lambda k: a.part in ('all', k)
    if run('hermite'):
        r = check_hermite()
        print("Hermite factorization of the eps-scores (four-cycle, m=6, s=1.3)")
        print("   tr H^2=%g tr H^3=%g tr H^4=%g  |H^3-4H|=%.0e |H^4-4H^2|=%.0e" % r["tr"])
        print("   L'''' rel=%.1e   f''/f rel=%.1e   f''''/f rel=%.1e\n"
              % (r["L4"], r["f2"], r["f4"]))
    if run('finite-diff'):
        print("-K_m/(m s^4) vs direct 4th differences of A_x")
        print("   %-4s %-5s %-6s %-17s %-17s %-9s"
              % ('m', 's', 'rho', 'theory', 'extrapolated FD', 'rel'))
        for row in check_against_finite_differences():
            print("   %-4d %-5s %-6s %-17.9g %-17.9g %-9.2e" % row)
        print()
    if run('table'):
        from curvature_conjectures import I_closed
        print("D^4 A at the two transitions")
        print("   %-4s %-5s %-12s %-15s %-12s %-15s"
              % ('m', 's', 'rho_c^-', 'D4A there', 'rho_c^+', 'D4A there'))
        for s in (1.0, 4.0):
            for m in (4, 6, 10):
                F = lambda r: I_closed(r, s, m-1)
                r1, r2 = _bisect(F, 0.34, 0.99), _bisect(F, 1.01, 12.0)
                print("   %-4d %-5s %-12.6f %-15.6g %-12.6f %-15.6g"
                      % (m, s, float(r1), float(D4A(s, float(r1)*s, m)),
                         float(r2), float(D4A(s, float(r2)*s, m))))
        print()
    if run('signs'):
        print("sign law: K_m(rho_c^-) > 0 and K_m(rho_c^+) < 0")
        print("   (the second half is FALSE at large m -- see the breakdown below)")
        rows, ok = sign_law()
        for m, s, r1, k1, r2, k2, good in rows:
            print("   m=%-4d s=%-6s rho_c^-=%-10.6f K=%-14.6g rho_c^+=%-10.6f "
                  "K=%-14.6g %s" % (m, s, r1, k1, r2, k2, "OK" if good else "FAIL"))
        print("   holds in every case:", ok)
        print()
        print("   K_m = (m-4)*Idd + 4*Icd at rho_c^+;  Idd = int W D^2 > 0 always")
        print("   %-4s %-5s %-11s %-14s %-14s %-14s"
              % ('m', 's', 'rho_c^+', 'Idd (>0)', 'Icd', 'K_m'))
        for row in sign_law_breakdown():
            print("   %-4d %-5s %-11.6f %-14.6g %-14.6g %-14.6g" % row)


if __name__ == '__main__':
    main()
