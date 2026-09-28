"""Feasible-query witnesses (plan Section 3).

A witness bounds min_q g_i(q) from ABOVE, never below. The value returned is therefore a
ceiling on the certified overlap any sound method can reach on these balls; a better
witness routine can only lower it (`witness_ceilings_bound_from_above_only`).
"""
import numpy as np


def overlap_adverse(score, member, k):
    """Members among the k smallest scores, ties resolved in favour of outsiders.

    Vectorised over leading axes. Equivalent to lexsort on (member, score) but O(n).
    """
    vk = np.partition(score, k - 1, axis=-1)[..., k - 1][..., None]
    lt = score < vk
    eq = score == vk
    n_lt = lt.sum(-1)
    mem_lt = (lt & member).sum(-1)
    out_eq = (eq & ~member).sum(-1)
    slots = k - n_lt
    return mem_lt + np.maximum(0, slots - out_eq)


def _sq_dists(Q, C):
    """Squared distances between query rows and centre rows, via BLAS."""
    d2 = (Q * Q).sum(1)[:, None] + (C * C).sum(1)[None, :] - 2.0 * (Q @ C.T)
    return np.maximum(d2, 0.0)


def _query_bank(c, rho, dc, E_idx, E_mask, i, n_targeted, n_random, rng):
    """Feasible queries in B_i: the centre, targeted directions, and boundary samples."""
    ci, ri, D = c[i], rho[i], c.shape[1]
    if ri == 0.0:
        return ci[None, :]
    dirs = [np.zeros(D)]

    di = dc[i].copy()
    di[i] = np.inf
    outs = np.flatnonzero(~E_mask[i])
    outs = outs[outs != i]
    near_out = outs[np.argsort(di[outs], kind="stable")[:n_targeted]]
    far_mem = E_idx[i][np.argsort(-di[E_idx[i]], kind="stable")[:n_targeted]]

    def unit(v):
        nv = np.linalg.norm(v)
        return v / nv if nv > 0 else v

    for j in near_out:                      # pull one outsider in
        dirs.append(unit(c[j] - ci))
    for l in far_mem:                       # push one member out
        dirs.append(-unit(c[l] - ci))
    if len(near_out) and len(far_mem):      # combined attack on the ranking
        v = sum(unit(c[j] - ci) for j in near_out) - sum(unit(c[l] - ci) for l in far_mem)
        dirs.append(unit(v))

    G = rng.normal(size=(n_random, D))
    G /= np.linalg.norm(G, axis=1, keepdims=True)
    G *= rng.uniform(0.6, 1.0, size=(n_random, 1))
    dirs.extend(list(G))
    return ci[None, :] + ri * np.asarray(dirs)


def witness_ceiling(c, rho, dc, E_idx, E_mask, k, rng,
                    n_targeted=8, n_random=32, batch=50):
    """v_i = min over evaluated feasible queries of g_i(q). Upper bound on H_i."""
    n, D = c.shape
    v = np.empty(n, dtype=int)
    max_slack = 0.0
    for lo in range(0, n, batch):
        hi = min(lo + batch, n)
        banks = [_query_bank(c, rho, dc, E_idx, E_mask, i, n_targeted, n_random, rng)
                 for i in range(lo, hi)]
        nq = max(len(b) for b in banks)
        Q = np.zeros((hi - lo, nq, D))
        valid = np.zeros((hi - lo, nq), bool)
        for a, b in enumerate(banks):
            Q[a, :len(b)] = b
            Q[a, len(b):] = b[0]                     # pad with the centre (feasible)
            valid[a, :len(b)] = True
        # verified ball membership (protocol: feasible_witness_membership_and_rank_verified)
        slack = np.linalg.norm(Q - c[lo:hi, None, :], axis=2) - rho[lo:hi, None]
        max_slack = max(max_slack, float(slack.max()))
        assert (slack <= 1e-9).all(), "infeasible witness query"

        d = np.sqrt(_sq_dists(Q.reshape(-1, D), c)).reshape(hi - lo, nq, n)
        mem = E_mask[lo:hi][:, None, :]
        score = np.where(mem, d + rho[None, None, :],
                         np.maximum(0.0, d - rho[None, None, :]))
        idx = np.arange(lo, hi)
        score[np.arange(hi - lo), :, idx] = np.inf   # exclude self
        g = overlap_adverse(score, np.broadcast_to(mem, score.shape), k)
        g = np.where(valid, g, k + 1)
        v[lo:hi] = g.min(axis=1)
    return v, max_slack
