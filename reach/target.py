"""The target of plan Section 3: reach, its box bounds, labels, and naive maps.

Ranks are 1-based positions in the stable order of distances from i, so ties are broken by
object ID as in the predecessor. Every distance matrix passed here has +inf on its diagonal.
"""
import numpy as np
from scipy.spatial.distance import cdist


def pairwise(A):
    """Euclidean distance matrix with +inf on the diagonal."""
    d = cdist(A, A)
    np.fill_diagonal(d, np.inf)
    return d


def rank_matrix(d):
    """rank[i, j] = 1-based position of j in the stable order of d[i]; self ranks last."""
    n = d.shape[0]
    order = np.argsort(d, axis=1, kind="stable")
    rank = np.empty((n, n), np.int32)
    np.put_along_axis(rank, order, np.arange(1, n + 1, dtype=np.int32)[None, :], axis=1)
    return rank


def m_of(eps, k):
    return int(np.ceil((1 - eps) * k - 1e-9))


def kth_smallest(A, m):
    """m-th smallest entry of every row (1-based m)."""
    return np.partition(A, m - 1, axis=1)[:, m - 1]


def reach_from_ranks(rank, E_idx, m):
    """reach_i = m-th smallest rank_i(j) over the displayed neighbors j in E_i."""
    return kth_smallest(np.take_along_axis(rank, E_idx, axis=1), m)


def displayed(emb, k):
    """E_idx (n x k) and E_mask (n x n) of a display's k nearest neighbors, stable ties."""
    de = pairwise(emb)
    n = len(emb)
    E_idx = np.argsort(de, axis=1, kind="stable")[:, :k]
    E_mask = np.zeros((n, n), bool)
    np.put_along_axis(E_mask, E_idx, True, axis=1)
    return E_idx, E_mask


def box_reach(L, U, E_idx, E_mask, m):
    """Exact reach extremes over the box L <= d <= U (plan Section 4).

    r_hi: members at U, outsiders at L; r_lo: the reverse. Valid with ID tie-breaking:
    an outsider ranked before the m-th member has d_il <= D_(m) <= U_(m), hence L_il <= U_(m);
    an outsider with U_il < L_(m) <= D_(m) is strictly closer, hence ranked before it.
    """
    n = L.shape[0]
    out = ~E_mask
    np.fill_diagonal(out, False)
    U_m = kth_smallest(np.take_along_axis(U, E_idx, axis=1), m)
    L_m = kth_smallest(np.take_along_axis(L, E_idx, axis=1), m)
    r_hi = m + ((L <= U_m[:, None]) & out).sum(1)
    r_lo = m + ((U < L_m[:, None]) & out).sum(1)
    return r_lo, r_hi


def member_reach(r_lo_j, r_hi_j, m):
    """Reach bounds from per-member rank bounds (comparison engines): m-th smallest of each."""
    return kth_smallest(r_lo_j, m), kth_smallest(r_hi_j, m)


def naive_ranks(c, s2, R, D, corrected):
    """Rank matrices of the naive maps: plug-in distances, or bias-corrected squared ones."""
    dc = pairwise(c)
    if not corrected:
        return rank_matrix(dc)
    d2 = dc ** 2 - (D / R) * s2[None, :]           # the row-constant -(D/R) s_i^2 drops out
    np.fill_diagonal(d2, np.inf)
    return rank_matrix(d2)
