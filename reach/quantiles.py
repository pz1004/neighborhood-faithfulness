"""Verified tail quantiles (plan Section 4).

Same policy as the predecessor's F quantile: the candidate comes from scipy, the guarantee
from evaluating the tail in arbitrary precision and moving the candidate outward until the
inequality is verified. Each returns a value q with the stated tail probability at most
alpha, so a bound built on it is conservative even if scipy's candidate is slightly off.
"""
import functools
import mpmath as mp
from scipy.stats import t as tdist, chi2

_STEP = mp.mpf("1e-12")


def _push(x, up):
    x = mp.mpf(x)
    return float(x * (1 + _STEP) + mp.mpf("1e-300")) if up else float(x * (1 - _STEP))


@functools.lru_cache(maxsize=None)
def t_upper(nu, alpha, dps=50):
    """Smallest verified t with P(T_nu > t) <= alpha (t > 0)."""
    t = float(tdist.isf(alpha, nu))
    with mp.workdps(dps):
        def sf(x):
            x = mp.mpf(x)
            return mp.betainc(mp.mpf(nu) / 2, mp.mpf(1) / 2, 0, nu / (nu + x * x),
                              regularized=True) / 2
        for _ in range(400):
            if sf(t) <= alpha:
                return t
            t = _push(t, up=True)
    raise RuntimeError("t quantile verification did not converge")


@functools.lru_cache(maxsize=None)
def chi2_upper(k, alpha, dps=50):
    """Smallest verified q with P(chi2_k > q) <= alpha."""
    q = float(chi2.isf(alpha, k))
    with mp.workdps(dps):
        def sf(x):
            return mp.gammainc(mp.mpf(k) / 2, mp.mpf(x) / 2, mp.inf, regularized=True)
        for _ in range(400):
            if sf(q) <= alpha:
                return q
            q = _push(q, up=True)
    raise RuntimeError("chi2 upper quantile verification did not converge")


@functools.lru_cache(maxsize=None)
def chi2_lower(k, alpha, dps=50):
    """Largest verified q with P(chi2_k <= q) <= alpha."""
    q = float(chi2.ppf(alpha, k))
    with mp.workdps(dps):
        def cdf(x):
            return mp.gammainc(mp.mpf(k) / 2, 0, mp.mpf(x) / 2, regularized=True)
        for _ in range(400):
            if cdf(q) <= alpha:
                return q
            q = _push(q, up=False)
    raise RuntimeError("chi2 lower quantile verification did not converge")
