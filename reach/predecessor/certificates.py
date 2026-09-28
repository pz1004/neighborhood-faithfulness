"""Certificate bounds of plan Section 4 that Stage 0 needs.

Implements the independent-distance baseline (4.1), the root-cell collective bound (4.2),
and the Section 4.5 clique-size ceiling with its sound early exit. No subdivision: Stage 0
runs without a certification solver (`studies.stage0_operating_region`).
"""
import numpy as np


def _topk_member_count(score, member, k):
    """Overlap of the displayed set with the k smallest scores, adverse ties.

    Ties are broken in favour of outsiders by making the member flag the secondary key.
    """
    order = np.lexsort((member.astype(np.int8), score), axis=-1)
    return np.take_along_axis(member, order[..., :k], axis=-1).sum(axis=-1)


def interval_bounds(dc, rho):
    """L_ij and U_ij of plan Section 4.1."""
    r = rho[:, None] + rho[None, :]
    return np.maximum(0.0, dc - r), dc + r


def h_box(dc, rho, E_mask, k):
    """Exact independent-interval overlap (plan Section 4.1); valid lower bound on H_i."""
    n = dc.shape[0]
    L, U = interval_bounds(dc, rho)
    score = np.where(E_mask, U, L)
    np.fill_diagonal(score, np.inf)
    return _topk_member_count(score, E_mask, k)


def clique_ceiling(t, rho_max, rho_min, cell_radius):
    """m_i(t,C) of plan Section 4.5: largest integer strictly below 1+2A^2/(B^2-2A^2).

    A and B use the sound uniform instantiation, so m is conservative (too large), which
    makes the early exit below conservative in the safe direction.
    """
    A = t + rho_max + cell_radius
    B = 2.0 * t + 2.0 * rho_min
    den = B * B - 2.0 * A * A
    if den <= 0:
        return np.inf
    x = 1.0 + 2.0 * A * A / den
    return max(1, int(np.ceil(x)) - 1)


def _greedy_clique_cover(adj, order):
    """Partition vertices into cliques of the conflict graph, greedily, in `order`.

    Every clique is validated: a vertex joins only if it conflicts with all current
    members (`certification.graph.validate_every_clique`). Unprocessed vertices become
    singletons, so the cover always covers every outsider.
    """
    cliques = []
    for v in order:
        for cl in cliques:
            if all(adj[v, u] for u in cl):
                cl.append(v)
                break
        else:
            cliques.append([v])
    return cliques


def collective_pass(dc, rho, E_idx, E_mask, k, build_graph=True):
    """One pass over all points: root-cell collective bound and Section 4.5 statistics.

    Returns h_graph (per point), a boolean `passes` (some threshold admits a nonzero
    collective bound), |O_i(t)| and m_i(t) at the binding threshold, and counters for
    graphs built versus early exits.
    """
    n = dc.shape[0]
    L, U = interval_bounds(dc, rho)
    hg = np.zeros(n, int)
    passes = np.zeros(n, bool)
    O_bind = np.zeros(n, int)
    m_bind = np.zeros(n, float)
    O_all = []
    built = exits = 0

    for i in range(n):
        out = np.flatnonzero(~E_mask[i])
        out = out[out != i]
        Lo = L[i][out]
        srt = np.argsort(Lo, kind="stable")
        Lo_s, out_s = Lo[srt], out[srt]
        rho_s = rho[out_s]
        # prefix extremes over the L-sorted outsiders: O_i(t) is always a prefix
        pmax = np.maximum.accumulate(rho_s) if len(rho_s) else np.zeros(0)
        pmin = np.minimum.accumulate(rho_s) if len(rho_s) else np.zeros(0)

        Ui = U[i][E_idx[i]]
        ts = np.unique(Ui)
        for j, t in enumerate(ts):
            s_t = int((Ui <= t).sum())
            nO = int(np.searchsorted(Lo_s, t, side="right"))
            if nO == 0:
                hg[i] = max(hg[i], min(s_t, k))
                passes[i] = True
                m = np.inf
            else:
                m = clique_ceiling(t, float(pmax[nO - 1]), float(pmin[nO - 1]), float(rho[i]))
                if np.isfinite(m) and nO >= k * m:
                    exits += 1                       # Section 4.5 early exit: h_graph(t)=0
                else:
                    passes[i] = True
                    if build_graph:
                        built += 1
                        O = out_s[:nO]
                        Rr = rho[O]
                        Lc = np.maximum(0.0, dc[np.ix_(O, O)] - Rr[:, None] - Rr[None, :])
                        adj = Lc > 2.0 * t
                        np.fill_diagonal(adj, False)
                        q = min(len(_greedy_clique_cover(adj, o)) for o in (
                            np.argsort(-adj.sum(1), kind="stable"),
                            np.arange(nO),
                            np.argsort(O, kind="stable")))
                        hg[i] = max(hg[i], min(s_t, max(0, k - q)))
            if j == len(ts) - 1:
                O_bind[i] = nO
                m_bind[i] = m
        O_all.append(O_bind[i])
    return hg, passes, O_bind, m_bind, built, exits
