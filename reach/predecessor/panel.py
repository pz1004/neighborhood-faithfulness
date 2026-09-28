"""One Stage-0 panel: generate, acquire, embed, and measure (plan Sections 5.4, 5.5)."""
import numpy as np
from scipy.spatial.distance import cdist

from . import generators as gen
from . import uncertainty as unc
from . import certificates as cert
from . import witness as wit
from .intrinsic import levina_bickel

EPSILONS = (0.10, 0.20, 0.30)


def pca2(c, dim=2):
    """One PCA fit per panel, on the observed means (plan Section 5.5)."""
    Z = c - c.mean(0)
    _, _, Vt = np.linalg.svd(Z, full_matrices=False)
    return Z @ Vt[:dim].T


def run_panel(family, gamma, n, k, R, delta, d0, seed):
    rng = np.random.default_rng(seed)
    mu = gen.draw_means(family, n, rng)

    exact_zero = (gamma == 0.0)
    X = unc.acquire(mu, gamma, d0, R, rng)
    c, rho, fq = unc.balls(X, delta, n, exact_zero=exact_zero)

    d_mu = cdist(mu, mu)
    np.fill_diagonal(d_mu, np.inf)
    N_star = np.argsort(d_mu, axis=1, kind="stable")[:, :k]

    emb = pca2(c)
    d_emb = cdist(emb, emb)
    np.fill_diagonal(d_emb, np.inf)
    E_idx = np.argsort(d_emb, axis=1, kind="stable")[:, :k]
    E_mask = np.zeros((n, n), bool)
    np.put_along_axis(E_mask, E_idx, True, axis=1)
    np.fill_diagonal(E_mask, False)

    realized_overlap = np.array([len(set(E_idx[i]) & set(N_star[i])) for i in range(n)])
    realized_err = 1.0 - realized_overlap / k

    dc = cdist(c, c)
    np.fill_diagonal(dc, np.inf)

    hb = cert.h_box(dc, rho, E_mask, k)
    hg, passes, O_bind, m_bind, built, exits = cert.collective_pass(
        dc, rho, E_idx, E_mask, k, build_graph=True)
    ell = np.maximum(hb, hg)

    v, max_slack = wit.witness_ceiling(c, rho, dc, E_idx, E_mask, k,
                                       np.random.default_rng(seed + 1))

    idim = levina_bickel(d_mu, k)

    assert (hb <= v).all(), "lower bound above witness ceiling"
    assert (hg <= v).all(), "collective bound above witness ceiling"
    if exact_zero:
        assert (hb == realized_overlap).all(), "zero-radius control: h_box must be exact"
        assert (v == realized_overlap).all(), "zero-radius control: ceiling must be exact"
        assert hg.max() == 0 or (hg <= hb).all(), "zero-radius control: no collective gain"

    row = dict(family=family, gamma=gamma, seed=seed, n=n, k=k, R=R, delta=delta,
               d0=d0, f_quantile=fq, rho_over_d0=float(np.mean(rho) / d0),
               witness_max_slack=max_slack,
               intrinsic_dim_median=float(np.nanmedian(idim)),
               realized_err_mean=float(realized_err.mean()),
               O_bind_median=float(np.median(O_bind)),
               O_bind_p90=float(np.percentile(O_bind, 90)),
               m_bind_median=float(np.median(m_bind[np.isfinite(m_bind)])
                                   if np.isfinite(m_bind).any() else np.inf),
               necessary_pass_rate=float(passes.mean()),
               graphs_built=built, early_exits=exits,
               graph_gain_points=int((hg > hb).sum()),
               graph_gain_mean=float(np.maximum(hg - hb, 0).mean()),
               ell_median=float(np.median(ell)), v_median=float(np.median(v)))
    for eps in EPSILONS:
        thr = int(np.ceil((1 - eps) * k))
        row[f"realized_ok@{eps}"] = float((realized_err <= eps + 1e-12).mean())
        row[f"cov_box@{eps}"] = float((hb >= thr).mean())
        row[f"cov_graph@{eps}"] = float((ell >= thr).mean())
        row[f"ceiling@{eps}"] = float((v >= thr).mean())
    return row, dict(ell=ell, hb=hb, hg=hg, v=v, realized=realized_overlap,
                     idim=idim, O_bind=O_bind, rho=rho)
