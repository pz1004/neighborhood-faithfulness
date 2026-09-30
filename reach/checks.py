"""Checks an analyst can run on a 2-D map from their own trials (docs/checks_protocol.md).

Each check gives every object a pass (True) or fail for its displayed neighborhood E_i.
"E_i recurs in a map M" means the reach of E_i computed from M's 2-D distances is at most
k': at least m of the k displayed neighbors are among i's k' nearest objects in M.
"""
import numpy as np
from sklearn.manifold import TSNE

from reach import target as tg

B_BOOT = 50                # bootstrap maps per display
BOOT_SUPPORT = 0.95        # BOOT passes if E_i recurs in at least this share of them
SCDEED_X = 0.5             # scDEED's default neighborhood, a share of all objects
SCDEED_PCT = 95.0          # scDEED's "trustworthy" cutoff, a percentile of the null
NAMES = ("N", "TW", "SCD", "CV", "SEED", "BOOT", "SH")


def tsne(c, seed, perplexity=30):
    """The primary cell's display: t-SNE with PCA initialization, one thread."""
    return TSNE(n_components=2, perplexity=perplexity, init="pca", random_state=seed,
                n_jobs=1).fit_transform(c)


def recurs(E_idx, emb, m, kp):
    """True where E_i recurs in the map `emb`."""
    return tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(emb)), E_idx, m) <= kp


def trustworthiness_points(rank, E_idx, k):
    """Per-point trustworthiness (Venna & Kaski); its mean over points is the usual score.

    T_i = 1 - 2 / (k (2n - 3k - 1)) * sum over j in E_i of max(0, rank_i(j) - k).
    """
    n = rank.shape[0]
    excess = np.maximum(np.take_along_axis(rank, E_idx, axis=1) - k, 0).sum(axis=1)
    return 1.0 - 2.0 / (k * (2 * n - 3 * k - 1)) * excess


def top_matched(score, n_pass):
    """Pass the n_pass highest scores, ties at the cutoff included; none if n_pass is 0."""
    if n_pass <= 0:
        return np.zeros(len(score), bool)
    cut = np.sort(score)[::-1][n_pass - 1]
    return score >= cut


def scdeed_scores(pre, emb, x=SCDEED_X):
    """scDEED reliability scores (Xia, Lee & Li 2024, Methods, step 4).

    For each object: the Pearson correlation of its 2-D distances to its round(x n) nearest
    objects in `pre`, in the order of their distance in `pre`, with its 2-D distances to its
    round(x n) nearest objects in the map, in the order of their 2-D distance.
    """
    n = len(pre)
    m = int(round(x * n))
    dpre, demb = tg.pairwise(pre), tg.pairwise(emb)
    nb = np.argsort(dpre, axis=1, kind="stable")[:, :m]
    v1 = np.take_along_axis(demb, nb, axis=1)
    v2 = np.sort(demb, axis=1)[:, :m]
    a = v1 - v1.mean(axis=1, keepdims=True)
    b = v2 - v2.mean(axis=1, keepdims=True)
    return (a * b).sum(axis=1) / np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1))


def scdeed_trustworthy(c, emb, seed, rng, c_raw, project, x=SCDEED_X, pct=SCDEED_PCT):
    """scDEED's "trustworthy" label (Methods, steps 1-6): score >= the pct-th null percentile.

    Null, as scDEED permutes the measured features before its PCA step: every raw feature of
    the averages `c_raw` is permuted independently across objects, the permuted averages go
    through the same fixed representation `project`, they are mapped with the display's
    t-SNE settings and seed, and all n of their scores form the null.
    """
    perm = project(np.column_stack([col[rng.permutation(len(col))] for col in c_raw.T]))
    null = scdeed_scores(perm, tsne(perm, seed), x)
    return scdeed_scores(c, emb, x) >= np.percentile(null, pct)


def bootstrap_support(Xu, E_idx, seed, m, kp, rng, B=None):
    """Share of B (default B_BOOT) bootstrap maps in which E_i recurs.

    Each resample draws every object's R trials with replacement, independently across
    objects, and is mapped with the display's t-SNE settings and seed.
    """
    B = B_BOOT if B is None else B
    n, R, _ = Xu.shape
    hits = np.zeros(n)
    for _ in range(B):
        idx = rng.integers(0, R, size=(n, R))
        cb = np.take_along_axis(Xu, idx[:, :, None], axis=1).mean(axis=1)
        hits += recurs(E_idx, tsne(cb, seed), m, kp)
    return hits / B


def labels_for_budget(Xu, c_raw, project, halves, seeds, rng_boot, rng_perm, cv_sq_dist,
                      k=10, eps=0.2):
    """Every check's labels for each display seed of one user budget.

    Xu: objects x R x D, the analyst's trials in the representation; c_raw: their averages
    in the measured features, and `project` the fixed map from those to the representation.
    `halves` are the two index sets of the RSA pipeline, `seeds` the display seeds, and
    rng_boot / rng_perm one generator per seed. Returns [(E_idx, {check: labels})].
    """
    m, kp = tg.m_of(eps, k), 2 * k
    A, B = halves
    c = Xu.mean(axis=1)
    rank = tg.rank_matrix(tg.pairwise(c))
    cv_rank = tg.rank_matrix(cv_sq_dist(Xu, A, B))
    cA, cB = Xu[:, A].mean(axis=1), Xu[:, B].mean(axis=1)
    embs = [tsne(c, s) for s in seeds]
    E = [tg.displayed(e, k)[0] for e in embs]
    out = []
    for j, s in enumerate(seeds):
        N = tg.reach_from_ranks(rank, E[j], m) <= kp
        lab = dict(N=N)
        lab["TW"] = top_matched(trustworthiness_points(rank, E[j], k), int(N.sum()))
        lab["SCD"] = scdeed_trustworthy(c, embs[j], s, rng_perm[j], c_raw, project)
        lab["CV"] = tg.reach_from_ranks(cv_rank, E[j], m) <= kp
        others = [embs[i] for i in range(len(seeds)) if i != j]
        lab["SEED"] = np.logical_and.reduce([recurs(E[j], e, m, kp) for e in others])
        lab["BOOT"] = bootstrap_support(Xu, E[j], s, m, kp, rng_boot[j]) >= BOOT_SUPPORT - 1e-9
        lab["SH"] = recurs(E[j], tsne(cA, s), m, kp) & recurs(E[j], tsne(cB, s), m, kp)
        out.append((E[j], lab))
    return out


def check_counts(lab, ref_reaches, kp):
    """Counts behind the joint reference rates of the protocol's Section 4, for one display.

    For every check X: X_pass, X_nx = #(N and X pass), and X_a = #(X passes, T_ref = 0) and
    X_b = #(N passes, X fails, T_ref = 0) summed over reference sub-blocks. For N, also the
    confirmatory counts n1, n10, dis and t1, so the rows can be compared with those rows.
    """
    N = lab["N"]
    row = dict(n=int(len(N)), nb=len(ref_reaches), n1=int(N.sum()), n10=0, dis=0, t1=0)
    for X in NAMES:
        row[f"{X}_pass"], row[f"{X}_nx"] = int(lab[X].sum()), int((N & lab[X]).sum())
        row[f"{X}_a"], row[f"{X}_b"] = 0, 0
    for rr in ref_reaches:
        T = rr <= kp
        row["n10"] += int((N & ~T).sum())
        row["dis"] += int((N != T).sum())
        row["t1"] += int(T.sum())
        for X in NAMES:
            S = lab[X]
            row[f"{X}_a"] += int((S & ~T).sum())
            row[f"{X}_b"] += int((N & ~S & ~T).sum())
    return row
