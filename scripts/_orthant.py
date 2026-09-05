"""Deterministic Genz separation-of-variables rectangle probabilities for
multivariate normals, with cancellation-free evaluation of upper tails.

P(max_j Y_j > x) is computed as the sum over "first index to exceed",
    sum_j P(Y_1<=x, ..., Y_{j-1}<=x, Y_j > x),
so no  1 - (1 - tiny)  subtraction ever occurs.  Nested Gauss-Legendre in
d-1 dimensions; deterministic, no Monte Carlo error.
"""
import numpy as np
from scipy.stats import norm

def _gl(N):
    t, w = np.polynomial.legendre.leggauss(N)
    return (t + 1)/2, w/2

def _seg(alpha, beta):
    """Phi(beta) - Phi(alpha), evaluated without cancellation."""
    both_pos = alpha > 0
    out = np.where(both_pos,
                   norm.sf(np.where(both_pos, alpha, 0.0)) - norm.sf(np.where(both_pos, beta, 0.0)),
                   norm.cdf(np.where(both_pos, 0.0, beta)) - norm.cdf(np.where(both_pos, 0.0, alpha)))
    return np.clip(out, 0.0, 1.0)

def rect_prob(a, b, L, N=64, max_nodes=4_000_000, N0=None):
    """P(a <= L z <= b), z ~ N(0, I_d), L lower-triangular with L[i,i] > 0.
    a, b may contain -inf / +inf."""
    d = len(a)
    if d > 1:
        N = max(16, min(N, int(round(max_nodes ** (1.0/(d-1))))))
    tw = {i: _gl(N) for i in range(d)}
    if N0: tw[0] = _gl(N0)          # the exceedance dimension needs many nodes
    def rec(i, y, wt, p):
        t, w = tw[i]
        off = (y @ L[i, :i]) if i else np.zeros(y.shape[0])
        alpha = np.atleast_1d((a[i] - off)/L[i, i])
        beta  = np.atleast_1d((b[i] - off)/L[i, i])
        seg = _seg(alpha, beta)
        p = p*seg
        if i == d - 1:
            return float(wt @ p)
        # invert in whichever tail is numerically safe: ppf(1-q) = -ppf(q)
        up = alpha > 0
        q_lo = np.where(up, norm.sf(np.where(up, alpha, 0.0)),
                            norm.cdf(np.where(up, 0.0, alpha)))
        grid = np.outer(seg, t)
        arg = np.clip(np.where(up[:, None], q_lo[:, None] - grid,
                                            q_lo[:, None] + grid), 1e-300, 1 - 1e-16)
        yi = np.where(up[:, None], -norm.ppf(arg), norm.ppf(arg))
        K, Ni = y.shape[0], len(t)
        return rec(i + 1,
                   np.concatenate([np.repeat(y, Ni, axis=0), yi.reshape(-1, 1)], axis=1),
                   np.repeat(wt, Ni)*np.tile(w, K), np.repeat(p, Ni))
    return rec(0, np.zeros((1, 0)), np.ones(1), np.ones(1))

def p_max_exceeds(mean, Sig, x, N=64, N0=600):
    """P(max_j Y_j > x) for Y ~ N(mean, Sig), computed tail-first."""
    d = len(mean)
    jit = 1e-14*max(1.0, float(np.max(np.abs(Sig))))
    u0 = np.asarray(x - mean, float)
    total = 0.0
    for j in range(d):                       # first index to exceed
        order = [j] + [k for k in range(j)]  # exceeding coord first, then Y_k<=x, k<j
        Sj = np.asarray(Sig, float)[np.ix_(order, order)]
        L = np.linalg.cholesky(Sj + jit*np.eye(len(order)))
        u = u0[order]
        a = np.full(len(order), -np.inf); b = u.copy()
        a[0] = u[0]; b[0] = np.inf
        total += rect_prob(a, b, L, N, N0=N0)
    return total

def p_max_exceeds_fast(mean, Sig, x, N=64):
    """P(max_j Y_j > x) as 1 - P(all <= x).  Fast, but loses absolute
    accuracy ~1e-10, so use only where the answer is >~ 1e-7."""
    d = len(mean)
    jit = 1e-14*max(1.0, float(np.max(np.abs(Sig))))
    L = np.linalg.cholesky(np.asarray(Sig, float) + jit*np.eye(d))
    return 1.0 - rect_prob(np.full(d, -np.inf), np.asarray(x - mean, float), L, N)
