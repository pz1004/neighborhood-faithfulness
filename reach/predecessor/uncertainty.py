"""Confidence balls of plan Section 2.1, with the verified F quantile of Section 9.1."""
import functools
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import f as fdist
import mpmath as mp


@functools.lru_cache(maxsize=None)
def verified_f_quantile(D, R, n, delta, dps=60):
    """Smallest verified f with  P(F_{D,D(R-1)} <= f) >= 1 - delta/n.

    Plan Section 9.1: certify by verification, not derivation. The candidate comes from
    scipy; the guarantee comes from evaluating the regularized incomplete beta in
    arbitrary precision and raising f until the inequality is verified. One scalar per
    configuration, computed once -- its cost enters no budget.
    """
    d1, d2 = D, D * (R - 1)
    target = 1 - mp.mpf(delta) / n
    f = float(fdist.ppf(1 - delta / n, d1, d2))
    with mp.workdps(dps):
        def cdf(x):
            x = mp.mpf(x)
            z = d1 * x / (d1 * x + d2)
            return mp.betainc(mp.mpf(d1) / 2, mp.mpf(d2) / 2, 0, z, regularized=True)
        for _ in range(200):
            if cdf(f) >= target:
                return float(f), True
            f = float(mp.mpf(f) * (1 + mp.mpf("1e-12")) + mp.mpf("1e-300"))
    raise RuntimeError("F quantile verification did not converge")


def acquire(mu, gamma, d0, R, rng):
    """Replicate acquisitions of plan Section 2.1; sigma_i = gamma*d0*w_i/sqrt(D), w_i=1."""
    n, D = mu.shape
    sigma = gamma * d0 / np.sqrt(D)
    X = mu[:, None, :] + sigma * rng.normal(size=(n, R, D))
    return X


def balls(X, delta, n_union, exact_zero=False):
    """Centres c_i and radii rho_i.

    exact_zero implements the zero-noise control: declared rho_i = 0, bypassing the
    positive-variance F pivot (`uncertainty.zero_noise_control`).
    """
    n, R, D = X.shape
    c = X.mean(axis=1)
    if exact_zero:
        return c, np.zeros(n), None
    resid = X - c[:, None, :]
    s2 = (resid ** 2).sum(axis=(1, 2)) / (D * (R - 1))
    q, verified = verified_f_quantile(D, R, n_union, delta)
    assert verified
    rho = np.sqrt(s2 * (D / R) * q)
    return c, rho, q


def reference_neighbors(mu, k):
    """N_i*: k nearest OTHER means, deterministic object-ID ties (plan Section 2)."""
    d = cdist(mu, mu)
    np.fill_diagonal(d, np.inf)
    return np.argsort(d, axis=1, kind="stable")[:, :k], d


def median_kth_neighbor_distance(mu, k):
    """d0 on a development panel (plan Section 7.1)."""
    d = cdist(mu, mu)
    np.fill_diagonal(d, np.inf)
    return float(np.median(np.sort(d, axis=1)[:, k - 1]))
