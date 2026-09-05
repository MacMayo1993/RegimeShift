"""Numerical evidence for the two open conjectures of the curvature-transition
paper (Conjecture 15: exactly two curvature zeros; Conjecture 16: global
scale-profiled cyclic maximality).

Nothing here is a proof. The script does four things, each reproducible and
free of Monte Carlo error:

1. ``check_theorem2`` -- confirms the closed three-integral form (1.11) against
   the raw Hessian integrand (1.8) to full working precision. This validates
   the algebra every other check below is built on.

2. ``scan_concavity`` -- evaluates the normalized curvature

       Q_{m,s}(rho) = I_m(s, rho*s) / (4 sqrt(s) J_1(b)),   b = sqrt(s)(1+rho),

   and its exact second derivative on an adversarial grid in (n, s, rho). The
   proposed sufficient condition for Conjecture 15 is Q''(rho) < 0 for every
   rho > 1/3; strict concavity there admits at most two zeros, and the three
   known facts Q(1/3) < 0, Q(1) > 0, Q -> -infinity force at least two.

3. ``check_identities`` -- verifies the two exact identities of the module
   docstring below, which are the intended route to a proof of (2).

4. ``profiled_family`` -- the m=4 one-parameter cyclic sweep for Conjecture 16,
   including the rank-deficient boundary profile w=1, which is handled by its
   own exact low-rank formula rather than by a Cholesky factor of a singular
   covariance.

Two exact identities are used and verified. Under the tilted density

    p_b(t) = phi(t-b) phi(t) Phi(t)^{n-1} / J_1(b),   R = phi/Phi,

    (A)   E_b[T] = b/2 + (n-1) u(b) / 2,          u = E_b[R],
    (B)   Cov_b(T, f(T)) = E_b[f'(T)]/2 + ((n-1)/2) Cov_b(f(T), R(T)).

(B) follows from integrating by parts against
(log p_b)'(t) = b - 2t + (n-1) R(t); (A) is (B) with f == 1. Together with
d/db E_b f = Cov_b(T, f) they express Q'' entirely in truncated-normal
moments, which is what ``Q_and_derivs`` computes.

Requires mpmath (for the curvature part) and numpy/scipy (for the profiled
part). Run ``python scripts/curvature_conjectures.py --help``.
"""

from mpmath import mp, mpf, exp, log, sqrt, quad, erfc, findroot

mp.dps = 30

def logphi(t):  return -t*t/2 - log(2*mp.pi)/2
def logPhi(t):  return log(erfc(-t/sqrt(2))/2)
def R(t):       return exp(logphi(t) - logPhi(t))

def _peak(b, j, n):
    f  = lambda t: -(t-b) - j*t + (n-j)*R(t)
    guess = (b + sqrt(2*log(max(n,2)))*(n-j)/mpf(n+1))/(1+j)
    try:    return findroot(f, mpf(guess))
    except Exception: return mpf(guess)

def _int(b, n, g, width=45):
    """integral of g(t) * phi(t-b) phi(t) Phi(t)^{n-1} dt"""
    t0 = _peak(b, 1, n)
    f = lambda t: g(t)*exp(-(t-b)**2/2 - t*t/2 + (n-1)*logPhi(t))/(2*mp.pi)
    return quad(f, [t0-width, t0-6, t0-1, t0+1, t0+6, t0+width])

def moments(b, n):
    """M[p][k] = E_b[ T^p R(T)^k ] under p_b, for p<=2, k<=2."""
    b = mpf(b)
    raw = {}
    for p in range(3):
        for k in range(3):
            raw[(p,k)] = _int(b, n, lambda t, p=p, k=k: t**p * R(t)**k)
    Z = raw[(0,0)]                      # = J_1(b)
    return {pk: val/Z for pk, val in raw.items()}, Z

def uv_derivs(b, n):
    """(u,u',u'', v,v',v'') in b, exact via Cov / second-central-moment identities."""
    M, _ = moments(b, n)
    m10, m20 = M[(1,0)], M[(2,0)]
    s2 = m20 - m10**2
    out = []
    for k in (1, 2):
        f0, f1, f2 = M[(0,k)], M[(1,k)], M[(2,k)]
        d1 = f1 - m10*f0                                  # Cov(T, f)
        d2 = (f2 - 2*m10*f1 + m10**2*f0) - s2*f0          # E[((T-mu)^2-sig^2) f]
        out += [f0, d1, d2]
    return out            # u,u',u'', v,v',v''

# ---- Q = A + B u + C v ----
def A(rho, s, n):  return (3*rho-1)/mpf(2) - s/mpf(4)*(1+rho)*(1-rho)**2
def B(rho, s, n):  return sqrt(s)/mpf(12)*(n*(3*rho**2+10*rho-5) - 3*(1-rho)**2)
def C(rho, s, n):  return n*(n+1)*(3*rho-1)/mpf(6)

def Q_and_derivs(rho, s, n):
    rho, s = mpf(rho), mpf(s)
    b = sqrt(s)*(1+rho)
    u, u1, u2, v, v1, v2 = uv_derivs(b, n)
    Ap  = mpf(3)/2 - s/mpf(4)*(3*rho**2 - 2*rho - 1)
    App = -(s/mpf(2))*(3*rho-1)
    Bp  = sqrt(s)/mpf(12)*(n*(6*rho+10) + 6*(1-rho))
    Bpp = sqrt(s)*(n-1)/mpf(2)
    Cp  = n*(n+1)/mpf(2)
    Bv, Cv = B(rho,s,n), C(rho,s,n)
    q0 = A(rho,s,n) + Bv*u + Cv*v
    q1 = Ap + Bp*u + sqrt(s)*Bv*u1 + Cp*v + sqrt(s)*Cv*v1
    q2 = (App + Bpp*u + 2*Bp*sqrt(s)*u1 + Bv*s*u2
              + 2*Cp*sqrt(s)*v1 + Cv*s*v2)
    parts = dict(App=App, Bpp_u=Bpp*u, cross_u=2*Bp*sqrt(s)*u1, curv_u=Bv*s*u2,
                 cross_v=2*Cp*sqrt(s)*v1, curv_v=Cv*s*v2)
    return q0, q1, q2, parts

def Q(rho, s, n):   return Q_and_derivs(rho, s, n)[0]
def Q2(rho, s, n):  return Q_and_derivs(rho, s, n)[2]

# ---- independent cross-check of Theorem 2 ----
def J(b, j, n, width=45):
    b = mpf(b); t0 = _peak(b, j, n)
    f = lambda t: exp(-(t-b)**2/2 - j*t*t/2 + (n-j)*logPhi(t))/(2*mp.pi)**((j+1)/mpf(2))
    return quad(f, [t0-width, t0-6, t0-1, t0+1, t0+6, t0+width])

def I_raw(rho, s, n, width=45):
    rho, s = mpf(rho), mpf(s); m = n+1; b = sqrt(s)*(1+rho)
    def integrand(t):
        r = R(t); Bt = -2*r*(t+r)
        At = (t-rho*sqrt(s))**2 - 1 + r*(t-2*rho*sqrt(s))
        return exp(-(t-b)**2/2 - log(2*mp.pi)/2 + n*logPhi(t))*((m-4)*Bt**2 + 4*At*Bt)
    t0 = _peak(b, 0, n)
    return quad(integrand, [t0-width, t0-6, t0-1, t0+1, t0+6, t0+width])

def I_closed(rho, s, n):
    rho, s = mpf(rho), mpf(s); b = sqrt(s)*(1+rho)
    c1 = sqrt(s)/mpf(4)*(2*(3*rho-1) - s*(1+rho)*(1-rho)**2)
    c2 = s/mpf(12)*(n*(3*rho**2+10*rho-5) - 3*(1-rho)**2)
    c3 = n*(n+1)*sqrt(s)*(3*rho-1)/mpf(6)
    return 4*(c1*J(b,1,n) + c2*J(b,2,n) + c3*J(b,3,n))


# ---------------------------------------------------------------- checks ----

def check_theorem2(cases=((0.5, 1, 3), (0.8, 4, 5), (1.5, 0.5, 9), (2.0, 9, 3))):
    """Closed form (1.11) against the raw Hessian integrand (1.8)."""
    print("Theorem 2: closed three-integral form vs raw integrand")
    ok = True
    for rho, s, n in cases:
        raw, closed = I_raw(rho, s, n), I_closed(rho, s, n)
        rel = abs(raw - closed)/abs(closed)
        ok &= rel < mpf(10)**(-20)
        print("  rho=%-5s s=%-5s n=%-4d  raw=%-22s closed=%-22s rel=%.1e"
              % (rho, s, n, mp.nstr(raw, 12), mp.nstr(closed, 12), float(rel)))
    return ok


def check_identities(cases=((0, 3), (1.3, 3), (2.7, 9), (-1.1, 5), (4, 100))):
    """Identities (A) and (B) of the module docstring."""
    print("Identities (A) E[T] = b/2 + (n-1)u/2   and   (B) the tilt identity")
    ok = True
    for b, n in cases:
        b = mpf(b)
        Z = _int(b, n, lambda t: mpf(1))
        E = lambda g: _int(b, n, g)/Z
        u, ET = E(R), E(lambda t: t)
        relA = abs(ET - (b/2 + (n-1)*u/2))/abs(ET)
        ok &= relA < mpf(10)**(-20)
        print("  n=%-4d b=%-6s  (A) rel=%.1e" % (n, mp.nstr(b, 4), float(relA)))
        for name, f, fp in [("R", R, lambda t: -R(t)*(t + R(t))),
                            ("R^2", lambda t: R(t)**2,
                             lambda t: -2*R(t)**2*(t + R(t)))]:
            cov = E(lambda t: t*f(t)) - ET*E(f)
            rhs = E(fp)/2 + (n-1)*(E(lambda t: f(t)*R(t)) - E(f)*u)/2
            rel = abs(cov - rhs)/abs(cov)
            ok &= rel < mpf(10)**(-20)
            print("      f=%-4s Cov(T,f)=%-20s rel=%.1e"
                  % (name, mp.nstr(cov, 12), float(rel)))
    return ok


DEFAULT_RHOS = ['1e-9', '1e-4', '0.01', '0.05', '0.15', '0.35',
                '0.667', '1.2', '2.0', '4.0', '8.0']
DEFAULT_NS = (3, 4, 5, 8, 20, 100, 1000, 9999)
DEFAULT_SS = ('1e-6', '1e-3', '0.03', '0.3', '1', '3', '10', '40')


def scan_concavity(ns=DEFAULT_NS, ss=DEFAULT_SS, rho_offsets=DEFAULT_RHOS,
                   verbose=True):
    """Search for a violation of the proposed sufficient condition Q'' < 0.

    Returns (rows, worst) where worst is the largest Q'' encountered. The grid
    deliberately includes rho within 1e-9 of 1/3 (where the polynomial part of
    Q'', which equals -(s/2)(3 rho - 1), vanishes and the sign is decided
    entirely by the integral terms), s down to 1e-6 and up to 40, and n up to
    9999.
    """
    rows, worst = [], None
    for n in ns:
        for s_str in ss:
            s = mpf(s_str)
            for off in rho_offsets:
                rho = mpf(1)/3 + mpf(off)
                q0, q1, q2, parts = Q_and_derivs(rho, s, n)
                rows.append(dict(n=n, s=s_str, rho=float(rho), Q=float(q0),
                                 Qp=float(q1), Qpp=float(q2),
                                 poly=float(parts['App']),
                                 positive_term=float(parts['Bpp_u'])))
                if worst is None or q2 > worst[0]:
                    worst = (q2, n, s_str, float(rho))
                if verbose:
                    print("  n=%-6d s=%-6s rho=%-10.5f Q=%-+14.6g Q''=%-+14.6g%s"
                          % (n, s_str, float(rho), float(q0), float(q2),
                             "   *** POSITIVE ***" if q2 > 0 else ""))
    return rows, worst


def small_s_limit(ns=(3, 4, 5, 10, 30, 100, 1000, 10000)):
    """As s -> 0 the two O(sqrt s) terms decide the sign of Q'':

        Q''(rho)/sqrt(s)  ->  L(n) = (n-1) u(0)/2 + n(n+1) v'(0),

    uniformly in rho, since the polynomial part is only O(s). The first term is
    positive -- it is the obstruction the user's note refers to -- and L(n) < 0
    only because n(n+1) v'(0) outweighs it.
    """
    out = []
    for n in ns:
        u, u1, u2, v, v1, v2 = uv_derivs(mpf(0), n)
        L = (n-1)*u/2 + n*(n+1)*v1
        out.append((n, u, v1, L))
        print("  n=%-6d u(0)=%-14s v'(0)=%-14s (n-1)u/2=%-12s L(n)=%-14s"
              % (n, mp.nstr(u, 6), mp.nstr(v1, 6), mp.nstr((n-1)*u/2, 6),
                 mp.nstr(L, 6)))
    return out


# -------------------------------------------- Conjecture 16: profiled tail ---
# numpy/scipy only from here down.

def _profiled_module():
    import numpy as np
    from math import cos, pi, sqrt as fsqrt, comb
    from scipy.stats import norm
    from scipy.integrate import quad as squad
    from scipy.optimize import minimize_scalar
    return np, cos, pi, fsqrt, comb, norm, squad, minimize_scalar


def sandwich(m, x):
    """The proved bracket

        n*Qb(sqrt(2x)) - C(n,2)*Qb(sqrt(8x/3)) <= F_x(G_tri) <= M_m(x) <= n*Qb(sqrt(2x))

    with Qb the standard normal survival function and n = m-1. The upper bound
    is a union bound valid over every unit-diagonal Gram matrix, not only the
    cyclic family; the lower bound is Bonferroni at the profiled simplex, where
    each competitor attains the individual tail exactly and each pair of
    standardized competitor scores has correlation 1/2.
    """
    np, cos, pi, fsqrt, comb, norm, squad, _ = _profiled_module()
    n = m - 1
    upper = n*norm.sf(fsqrt(2*x))
    return upper - comb(n, 2)*norm.sf(fsqrt(8*x/3)), upper


def gram(w, m):
    """Circulant Gram matrix of the centered cyclic profile w (length m//2)."""
    np, cos, pi, *_ = _profiled_module()
    q = m // 2
    return np.array([[sum(w[k-1]*cos(2*pi*k*(bb-a)/m) for k in range(1, q+1))
                      for bb in range(m)] for a in range(m)])


def Ax_profile(G, E, x, m, N=64):
    """A_x(EG); circulant G makes every transmitted index equivalent.

    Uses 1 - P(all competitors <= x), which is accurate only while the answer
    is above roughly 1e-7. Deep-tail values need the low-rank/tail-first route.
    """
    np, *_ = _profiled_module()
    from _orthant import p_max_exceeds_fast
    idx = list(range(1, m))
    mean = np.array([-E*(1 - G[0, j]) for j in idx])
    Sig = np.array([[E*(G[i, j] - G[0, i] - G[0, j] + 1) for j in idx]
                    for i in idx])
    keep = np.diag(Sig) > 1e-13*max(1.0, E)
    if not keep.any():
        return 0.0
    return p_max_exceeds_fast(mean[keep], Sig[np.ix_(keep, keep)], x, N)


def profiled(G, x, m, N=64):
    """F_x(G) = sup_{E>0} A_x(EG), by multistart bounded search in log E."""
    np, *rest = _profiled_module()
    minimize_scalar = rest[-1]
    f = lambda lE: -Ax_profile(G, np.exp(lE), x, m, N)
    best = min((minimize_scalar(f, bounds=bb, method='bounded',
                                options=dict(xatol=1e-9))
                for bb in [(-6, -1), (-1, 1.5), (1, 3), (2.5, 6)]),
               key=lambda r: r.fun)
    return -best.fun, float(np.exp(best.x))


def F_biorthogonal(x):
    """F_x at the m=4 boundary profile w=1, computed exactly.

    That Gram has rank 2 -- the signals are (1,0),(0,1),(-1,0),(0,-1) -- so the
    competitor covariance is singular and no Cholesky-based orthant routine
    applies. Writing W, V for iid N(0,E) coordinates of the noise,

        max_j l_0j > x   iff   W + |V| > x + E   or   W > x/2 + E,

    which is a one-dimensional quadrature. Ignoring the degeneracy and
    regularizing the singular covariance instead produces a spurious upward
    jump at w=1 that looks like a second local maximum.
    """
    np, cos, pi, fsqrt, comb, norm, squad, minimize_scalar = _profiled_module()

    def A(E):
        sd = fsqrt(E)
        a, c = x + E, x/2 + E
        f = lambda w: norm.pdf(w, 0, sd)*(2*norm.sf(a - w, 0, sd)
                                          if a - w > 0 else 1.0)
        return norm.sf(c, 0, sd) + squad(f, -12*sd, c, limit=400)[0]

    r = minimize_scalar(lambda lE: -A(np.exp(lE)), bounds=(-2, 3.5),
                        method='bounded', options=dict(xatol=1e-10))
    return -r.fun, float(np.exp(r.x))


def wstar(m):
    np, *_ = _profiled_module()
    n, q = m - 1, m // 2
    w = np.full(q, 2.0/n)
    if m % 2 == 0:
        w[-1] = 1.0/n
    return w


def profiled_family(m=4, xs=(0.5, 1.0, 2.0, 4.0), step=0.05):
    """Sweep the one-parameter m=4 cyclic family against the simplex."""
    np, *_ = _profiled_module()
    results = []
    for x in xs:
        Fstar, Estar = profiled(gram(wstar(m), m), x, m)
        lo, up = sandwich(m, x)
        grid = []
        for w in np.round(np.arange(0.0, 1.0, step), 6):
            F, E = profiled(gram([w, 1-w], m), x, m)
            grid.append((float(w), F, E))
        Fb, Eb = F_biorthogonal(x)
        grid.append((1.0, Fb, Eb))
        best = max(grid, key=lambda t: t[1])
        results.append(dict(x=x, F_simplex=Fstar, E_simplex=Estar,
                            E_theory=x*(m-1)/m, best_w=best[0], best_F=best[1],
                            bonferroni_lower=lo, union_upper=up, grid=grid))
        print("  x=%-5s F(w*)=%-14.9g E*=%-10.6g (theory %-10.6g) "
              "best off-simplex=%-14.9g at w=%-6s  bracket=[%.6g, %.6g] %s"
              % (x, Fstar, Estar, x*(m-1)/m,
                 max(g[1] for g in grid if abs(g[0] - 2/3) > 1e-9), best[0],
                 lo, up, "OK" if lo - 1e-12 <= Fstar <= up + 1e-12 else "OUT"))
    return results


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--part', default='all',
                   choices=['all', 'theorem2', 'identities', 'concavity',
                            'small-s', 'profiled'])
    p.add_argument('--quick', action='store_true',
                   help='small grids; the full concavity scan takes ~20 min')
    a = p.parse_args(argv)
    run = lambda k: a.part in ('all', k)
    if run('theorem2'):
        check_theorem2(); print()
    if run('identities'):
        check_identities(); print()
    if run('small-s'):
        print("Small-s limit  Q''/sqrt(s) -> L(n):")
        small_s_limit((3, 10, 100) if a.quick else None
                      or (3, 4, 5, 10, 30, 100, 1000, 10000)); print()
    if run('concavity'):
        print("Concavity scan: looking for any rho > 1/3 with Q'' >= 0")
        ns = (3, 20) if a.quick else DEFAULT_NS
        ss = ('1e-3', '1', '10') if a.quick else DEFAULT_SS
        rows, worst = scan_concavity(ns, ss)
        print("  max Q'' over %d points: %.6g at n=%s s=%s rho=%s"
              % (len(rows), float(worst[0]), worst[1], worst[2], worst[3]))
        print("  violations: %d" % sum(r['Qpp'] >= 0 for r in rows)); print()
    if run('profiled'):
        print("Conjecture 16, m=4 cyclic family:")
        profiled_family(4, (0.5, 2.0) if a.quick else (0.5, 1.0, 2.0, 4.0),
                        0.1 if a.quick else 0.05)


if __name__ == '__main__':
    main()
