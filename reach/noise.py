"""Noise models of plan Section 5.2. Every model has unit variance per coordinate at scale 1."""
import numpy as np

MODELS = ("gauss", "gauss_het", "t5", "t3", "rep_t5", "aniso4", "aniso16")
ANISO_SEED = 3_400_000


def _aniso_factor(kappa, D):
    """Sigma_0^{1/2} = Q diag(sqrt(lambda)): lambda geometric, max/min = kappa, mean 1."""
    lam = kappa ** np.linspace(0.0, 1.0, D)
    lam /= lam.mean()
    Q, _ = np.linalg.qr(np.random.default_rng(ANISO_SEED + kappa).normal(size=(D, D)))
    return Q * np.sqrt(lam)[None, :]          # columns scaled: e = z @ (Q diag(sqrt lam)).T


def acquire(mu, R, sigma, model, rng):
    """R replicates of every object: X[i, r] = mu_i + e_ir, with e from `model` at scale sigma."""
    n, D = mu.shape
    shape = (n, R, D)
    if model == "gauss":
        e = rng.normal(size=shape)
    elif model == "gauss_het":
        w = np.exp(rng.uniform(-np.log(2.0), np.log(2.0), n))
        e = rng.normal(size=shape) * w[:, None, None]
    elif model in ("t5", "t3"):
        nu = int(model[1])
        e = rng.standard_t(nu, size=shape) / np.sqrt(nu / (nu - 2.0))
    elif model == "rep_t5":
        nu = 5
        g = (nu - 2.0) / rng.chisquare(nu, size=(n, R, 1))
        e = rng.normal(size=shape) * np.sqrt(g)
    elif model in ("aniso4", "aniso16"):
        A = _aniso_factor(int(model[5:]), D)
        e = rng.normal(size=shape) @ A.T
    else:
        raise KeyError(model)
    return mu[:, None, :] + sigma * e
