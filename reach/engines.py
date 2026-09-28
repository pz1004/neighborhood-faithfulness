"""Engines of plan Section 4.

Box engines return intervals [L, U] for every pairwise distance, jointly valid with
probability at least 1 - delta. The comparison engine returns, for chosen members j of each
row, bounds on the true rank of j derived from certified pairwise comparisons.
"""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import norm, chi2

from .predecessor import uncertainty as unc
from .quantiles import t_upper, chi2_upper, chi2_lower

DELTA = 0.05
SPLIT = dict(dZ=0.040, dW=0.005, dS=0.005)       # plan Section 4, E-pair-t


def moments(X):
    """Replicate means c, pooled-over-coordinates variances s2, and nu = D(R-1)."""
    n, R, D = X.shape
    c = X.mean(axis=1)
    s2 = ((X - c[:, None, :]) ** 2).sum(axis=(1, 2)) / (D * (R - 1))
    return c, s2, D * (R - 1)


def ball(X, delta=DELTA):
    """E-ball: the predecessor's F-pivot balls turned into distance intervals."""
    c, rho, q = unc.balls(X, delta, X.shape[0])
    dc = cdist(c, c)
    r = rho[:, None] + rho[None, :]
    return dict(L=np.maximum(0.0, dc - r), U=dc + r, c=c, rho=rho, q=q)


def pair(X, pooled=False, dZ=SPLIT["dZ"], dW=SPLIT["dW"], dS=SPLIT["dS"]):
    """E-pair-t (pooled=False) or E-pair-pool (pooled=True).

    U_ij = |c_i - c_j| + t tau_hat_ij
    L_ij = max(0, sqrt(max(0, |c_i - c_j|^2 - q_W tau_bar_ij^2)) - t tau_hat_ij)
    tau_hat^2 = (s_i^2 + s_j^2)/R studentizes the projected noise (Student quantile,
    conservative for any sigma_i, sigma_j by the Jensen argument of plan Section 4);
    tau_bar^2 = (b_i + b_j)/R uses upper bounds b on sigma^2 for the orthogonal part.
    """
    n, R, D = X.shape
    c, s2, nu = moments(X)
    pairs = n * (n - 1) // 2
    if pooled:
        s2 = np.full(n, s2.mean())
        nu = n * nu
        b = s2 * nu / chi2_lower(nu, dS)                   # one event: the pooled variance
    else:
        b = s2 * nu / chi2_lower(nu, dS / n)               # n events, one per object
    t = t_upper(nu, dZ / (2 * pairs))
    qW = chi2_upper(D - 1, dW / pairs)
    dc = cdist(c, c)
    tau_hat = np.sqrt((s2[:, None] + s2[None, :]) / R)
    tau_bar2 = (b[:, None] + b[None, :]) / R
    U = dc + t * tau_hat
    L = np.maximum(0.0, np.sqrt(np.maximum(0.0, dc ** 2 - qW * tau_bar2)) - t * tau_hat)
    return dict(L=L, U=U, c=c, s2=s2, t=t, qW=qW, b=b, half_width=float(np.median(t * tau_hat)))


def cmp_member_ranks(X, J, delta=DELTA, mu=None):
    """E-cmp, the pilot's conservative comparison engine (reference only; approximate).

    For query i, member j = J[i, a] and any other l, the bias-corrected statistic
        Delta_ijl = |c_i-c_l|^2 - |c_i-c_j|^2 - (D/R)(s_l^2 - s_j^2)
    is compared with z * se, se from sigma^2 upper bounds and ball-inflated distances, z a
    Gaussian quantile Bonferroni-corrected over all n*C(n-1,2) comparisons (delta split
    three ways). Returns r_lo, r_hi (n x |J|): 1 + #{l certified closer} and
    1 + #{l not certified farther}. With `mu`, also counts false elementary statements.
    """
    n, R, Dm = X.shape
    c, s2, df = moments(X)
    dc = cdist(c, c)
    d2 = dc ** 2 - (Dm / R) * (s2[:, None] + s2[None, :])
    s2v = s2 * df / chi2.ppf(delta / 3 / n, df)
    _, rho, _ = unc.balls(X, delta / 3, n)
    dv2 = (dc + rho[:, None] + rho[None, :]) ** 2
    z = norm.isf(delta / 3 / (2 * n * (n - 1) * (n - 2) / 2))
    dmu2 = cdist(mu, mu) ** 2 if mu is not None else None
    kJ = J.shape[1]
    r_lo = np.empty((n, kJ), np.int32)
    r_hi = np.empty((n, kJ), np.int32)
    false_stmt = 0
    ar = np.arange(kJ)
    for i in range(n):
        Ji = J[i]
        Dl = d2[i][None, :] - d2[i, Ji][:, None]               # > 0: l farther than j
        var = 4.0 / R * (s2v[i] * dv2[Ji, :] + s2v[None, :] * dv2[i][None, :]
                         + s2v[Ji][:, None] * dv2[i, Ji][:, None])
        var += (2 * Dm / R ** 2) * (1 + 1 / (R - 1)) * (s2v[None, :] ** 2 + s2v[Ji][:, None] ** 2)
        var += 4 * Dm * s2v[i] * (s2v[None, :] + s2v[Ji][:, None]) / R ** 2
        se = np.sqrt(var)
        ok = np.ones((kJ, n), bool)
        ok[:, i] = False
        ok[ar, Ji] = False
        farther = ((Dl - z * se) > 0) & ok
        closer = ((Dl + z * se) < 0) & ok
        r_hi[i] = 1 + (ok & ~farther).sum(1)
        r_lo[i] = 1 + closer.sum(1)
        if dmu2 is not None:
            T = dmu2[i][None, :] - dmu2[i, Ji][:, None]
            false_stmt += int((farther & (T <= 0)).sum() + (closer & (T >= 0)).sum())
    return r_lo, r_hi, dict(z=float(z), false_statements=false_stmt)
