"""Simulation study of the false-comfort estimator (docs/revision_protocol.md, Section 3).

Every mechanism is built in the measured features of one real subject's split-0 blocks, its
*base*: the truth families T1-T3, noise resampled from block U's within-image residuals at a
scale per block, shared offsets per image and session or stratum (mechanism D), and an
analyst's block drawn with less noise than the rest (the stress test S).

One repetition is one virtual subject:
- a new dataset of the real design's size (80 trials per image in EEG, 50 in mouse V1);
- the published procedure on it, run with new splits:
  - the representation fitted on block F of each split;
  - the maps, with the published display seeds;
  - the naive labels and the references;
  - exact labels against the truth, projected through each split's representation.

Counts are pooled within the repetition, as the analyses pool within a subject. Layer A
regenerates only the references of the first split's maps.
"""
import numpy as np

from reach import things, allen, target as tg
import run_things_confirmatory as RT
import run_allen_confirmatory as RA

K, EPS = 10, 0.2
M, KP = tg.m_of(EPS, K), 2 * K
LAYOUT = {
    "EEG": dict(RU=(8, 16, 32), RH=(8, 20, 40), H=40, extra=10, deep=(80, 160), n_strata=4,
                strata=things.sessions, split=things.split_blocks, halves=things.halves,
                split_base=RT.SPLIT_BASE),
    "mouse": dict(RU=(8, 16), RH=(8, 12, 24), H=24, extra=6, deep=(48, 96), n_strata=2,
                  strata=allen.strata, split=allen.split_blocks, halves=allen.halves,
                  split_base=RA.SPLIT_BASE),
}
FIELDS = ("n1nb", "n10", "dis", "t1", "nnb", "n1", "n")     # pooled counts stored per reference


def load_feat(ds, s):
    """Subject s's responses in the primary window (images x repetitions x features), digest-checked."""
    if ds == "EEG":
        X, _, times = RT.verified_load(s)
        return things.features_window(X, times, *RT.WINDOWS["full"]).astype(np.float64)
    feats, _ = RA.verified_load(RA.included()[s - 1])
    return feats["full"]


class Base:
    """One real subject's split-0 blocks: the T1 truth, the residual pool and the layout."""

    def __init__(self, ds, s, feat=None):
        self.ds, self.s, self.L = ds, s, LAYOUT[ds]
        feat = load_feat(ds, s) if feat is None else feat
        self.strata = self.L["strata"]()
        b0 = self.L["split"](self.strata, np.random.default_rng(self.L["split_base"] + 1000 * s))
        U = feat[:, b0["U"]]
        nU = U.shape[1]
        self.mu1 = feat[:, b0["H"]].mean(axis=1)
        self.pool = (U - U.mean(axis=1, keepdims=True)) * np.sqrt(nU / (nU - 1))
        self.W0, self.m00 = RT.fit_variant("white", feat[:, b0["F"]])     # for the gap of T3
        self.n, self.reps, self.p = feat.shape


def truth(base, family, r=5, alpha=None):
    """The truth families of protocol Section 3.2, in the measured features."""
    mu1 = base.mu1
    if family == "T1":
        return mu1
    m = mu1.mean(axis=0)
    V = np.linalg.svd(mu1 - m, full_matrices=False)[2][:r].T
    mu2 = m + (mu1 - m) @ V @ V.T
    if family == "T2":
        return mu2
    assert family == "T3" and alpha is not None
    return mu2 + alpha * (mu1 - mu2)


def relative_gap(base, mu):
    """Median over images of (d_(20) - d_(8)) / d_(8) of the true distances, in the base split's
    representation: the quantity explore_*_snr.py estimates from cross-validated distances."""
    srt = np.sort(tg.pairwise((mu - base.m00) @ base.W0), axis=1)     # the diagonal's +inf sorts last
    return float(np.median((srt[:, KP - 1] - srt[:, M - 1]) / srt[:, M - 1]))


def draw_noise(base, rng, n_trials):
    """Unit-scale noise of n_trials trials of every image, resampled from the image's pool."""
    idx = rng.integers(0, base.pool.shape[1], size=(base.n, n_trials))
    return np.take_along_axis(base.pool, idx[:, :, None], axis=1)


def mean_noise(base, rng, n_trials):
    """The mean of n_trials resampled residuals per image, drawn without forming the trials."""
    nU = base.pool.shape[1]
    cnt = rng.multinomial(n_trials, np.full(nU, 1.0 / nU), size=base.n)
    return np.einsum("iu,iup->ip", cnt, base.pool) / n_trials


def halves_of(size, n_strata):
    """Two halves of a block in round-robin order: alternate rows where that balances the strata,
    otherwise even and odd positions."""
    p = np.arange(size)
    if size % (2 * n_strata) == 0:
        return p[(p // n_strata) % 2 == 0], p[(p // n_strata) % 2 == 1]
    return p[0::2], p[1::2]


def reference_ranks(YH, L):
    """Rank matrices of the references from block H (in round-robin order): plug-in distances of
    sub-block means, and cross-validated distances between the two halves of each sub-block, at
    the published reference sizes and at the extra size of candidate (b)."""
    out = {}
    for Rh in L["RH"] + (L["extra"],):
        subs = [YH[:, j * Rh:(j + 1) * Rh] for j in range(L["H"] // Rh)]
        out[("plug", Rh)] = [tg.rank_matrix(tg.pairwise(S.mean(axis=1))) for S in subs]
        A, B = halves_of(Rh, L["n_strata"])
        out[("cv", Rh)] = [tg.rank_matrix(things.cv_sq_dist(S, A, B)) for S in subs]
    return out


class Repetition:
    """One virtual subject: a dataset drawn from a cell's mechanism, and the published procedure on it.

    cell: dict(name, family, map ('tsne' or 'mds_cv'), R_list, lam_F, lam_U, lam_H, rho).
    """

    def __init__(self, base, cell, mu, rng_data):
        self.base, self.cell, self.mu = base, cell, mu
        rho = cell.get("rho", 0.0)
        E = draw_noise(base, rng_data, base.reps)
        self.O = None
        if rho > 0:
            self.O = draw_noise(base, rng_data, base.L["n_strata"])          # one residual per image and stratum
            E = np.sqrt(1 - rho) * E + np.sqrt(rho) * self.O[:, base.strata]
        self.E = E
        # deep references, beyond the design: their means only, with the same mechanism as block H
        self.deep_raw = {Rh: mu + cell["lam_H"] * self._mean_noise(rng_data, Rh) for Rh in base.L["deep"]}

    def _mean_noise(self, rng, n_trials):
        """Mean noise of n_trials trials spread evenly over the strata, with the shared offsets."""
        rho = self.cell.get("rho", 0.0)
        e = mean_noise(self.base, rng, n_trials)
        return e if rho == 0 else np.sqrt(1 - rho) * e + np.sqrt(rho) * self.O.mean(axis=1)

    def split(self, rng_split):
        """Blocks of one split, the dataset with each block's noise scale, and the representation."""
        c, L = self.cell, self.base.L
        b = L["split"](self.base.strata, rng_split)
        scale = np.empty(self.base.reps)
        scale[b["F"]], scale[b["U"]], scale[b["H"]] = c["lam_F"], c["lam_U"], c["lam_H"]
        X = self.mu[:, None, :] + scale[None, :, None] * self.E
        W, m0 = RT.fit_variant("white", X[:, b["F"]])
        return b, X, W, m0

    def maps(self, YU, R):
        """The displays of the analyst's first R trials, with their naive labels: [(E_idx, N)]."""
        L, Xu = self.base.L, YU[:, :R]
        c = Xu.mean(axis=1)
        if self.cell["map"] == "tsne":
            nrank = tg.rank_matrix(tg.pairwise(c))
            embs = [RT.embed("tsne30", c, sd) for sd in RT.TSNE_SEEDS]
        else:
            d2cv = things.cv_sq_dist(Xu, *L["halves"](R))
            nrank = tg.rank_matrix(d2cv)
            embs = [RT.embed("mds_cv", c, sd, d2cv) for sd in RT.MDS_SEEDS]
        out = []
        for emb in embs:
            E_idx, _ = tg.displayed(emb, K)
            out.append((E_idx, tg.reach_from_ranks(nrank, E_idx, M) <= KP))
        return out


def add_counts(acc, key, N, ranks, E_idx):
    """Add one display's counts against a list of reference rank matrices to acc[key]."""
    c = RT.counts(N, [tg.reach_from_ranks(rk, E_idx, M) for rk in ranks], KP)
    a = acc.setdefault(key, np.zeros(len(FIELDS)))
    a += [c["n1"] * c["nb"], c["n10"], c["dis"], c["t1"], c["n"] * c["nb"], c["n1"], c["n"]]


def run_repetition(base, cell, mu, seed, n_splits=10, n_regen=20):
    """Layers B and A for one virtual subject. seed is a tuple; every draw derives from it.

    Returns rows: per budget R, Layer B's counts pooled over splits and seeds for every reference
    ('plug', 'cv', 'deep', 'exact'), and per regeneration of the first split's references, Layer A's.
    """
    rep = Repetition(base, cell, mu, np.random.default_rng(seed + (0,)))
    accB = {}
    first = None
    for j in range(n_splits):
        b, X, W, m0 = rep.split(np.random.default_rng(seed + (1, j)))
        proj = lambda A, m0=m0, W=W: (A - m0) @ W            # bound now: Layer A uses split 0's
        Y = proj(X)
        YU, YH = Y[:, b["U"]], Y[:, b["H"]]
        refs = reference_ranks(YH, base.L)
        refs.update({("deep", Rh): [tg.rank_matrix(tg.pairwise(proj(m)))] for Rh, m in rep.deep_raw.items()})
        refs[("exact", 0)] = [tg.rank_matrix(tg.pairwise(proj(mu)))]
        for R in cell["R_list"]:
            mps = rep.maps(YU, R)
            for E_idx, N in mps:
                for key, ranks in refs.items():
                    add_counts(accB, (R,) + key, N, ranks, E_idx)
            if j == 0 and R == cell["R_list"][0]:
                first = (b, proj, {R: mps}, refs[("exact", 0)])
    rows = [dict(layer="B", regen=-1, R=R, ref=ref, Rh=Rh, maps=n_splits * 3, **dict(zip(FIELDS, a.tolist())))
            for (R, ref, Rh), a in sorted(accB.items())]
    if n_regen:
        rows += layer_a(rep, first, np.random.default_rng(seed + (2,)), n_regen)
    return rows


def layer_a(rep, first, rng, n_regen):
    """Layer A: the first split's maps at the first budget, fixed; block H and the deep references
    regenerated n_regen times from the same mechanism."""
    base, c = rep.base, rep.cell
    b, proj, maps, exact_rank = first
    (R, mp), = maps.items()
    st_H = base.strata[b["H"]]
    rows = []
    acc = {}
    for E_idx, N in mp:
        add_counts(acc, (R, "exact", 0), N, exact_rank, E_idx)
    rows += [dict(layer="A", regen=-1, R=k[0], ref=k[1], Rh=k[2], maps=len(mp), **dict(zip(FIELDS, a.tolist())))
             for k, a in acc.items()]
    for g in range(n_regen):
        E_H = draw_noise(base, rng, base.L["H"])
        if c.get("rho", 0.0) > 0:
            E_H = np.sqrt(1 - c["rho"]) * E_H + np.sqrt(c["rho"]) * rep.O[:, st_H]
        refs = reference_ranks(proj(rep.mu[:, None, :] + c["lam_H"] * E_H), base.L)
        for Rh in base.L["deep"]:
            refs[("deep", Rh)] = [tg.rank_matrix(tg.pairwise(proj(rep.mu + c["lam_H"] * rep._mean_noise(rng, Rh))))]
        acc = {}
        for E_idx, N in mp:
            for key, ranks in refs.items():
                add_counts(acc, (R,) + key, N, ranks, E_idx)
        rows += [dict(layer="A", regen=g, R=k[0], ref=k[1], Rh=k[2], maps=len(mp), **dict(zip(FIELDS, a.tolist())))
                 for k, a in sorted(acc.items())]
    return rows
