"""Levina-Bickel local intrinsic dimension (plan H5; lit_survey reference [17])."""
import numpy as np


def levina_bickel(dist, k):
    """Per-point MLE m_k(x_i) = [ (1/(k-1)) sum_{j<k} log(T_k/T_j) ]^{-1}.

    `dist` is the full pairwise distance matrix on the reference means; the diagonal is
    ignored. Points whose k-th neighbour is at distance 0 (exact duplicates) return nan.
    """
    d = dist.copy()
    np.fill_diagonal(d, np.inf)
    T = np.sort(d, axis=1)[:, :k]                       # T[:,0]..T[:,k-1]
    Tk = T[:, k - 1][:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.log(Tk / T[:, : k - 1])
        inv = logs.sum(axis=1) / (k - 1)
        m = np.where(inv > 0, 1.0 / inv, np.nan)
    return m
