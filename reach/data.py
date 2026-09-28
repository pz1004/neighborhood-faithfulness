"""Development geometries of plan Section 5.1. Only data bundled with installed libraries."""
import numpy as np
from sklearn.datasets import load_digits

from .predecessor import generators as gen
from .predecessor.uncertainty import median_kth_neighbor_distance

K, D = 20, 20
N_DIGITS = 1000
DIGITS_SPLIT_SEED = 20260923            # the predecessor's split, reused to reproduce its panel
SURFACE_DEV_SEED, SURFACE_LAYOUT_SEED = 3_300_000, 3_310_000
SYN_DEV_BASE, SYN_LAYOUT_BASE = 3_300_100, 3_310_100     # + index in gen.PRIMARY (power runs)


def digits_panel():
    """The predecessor's real-data panel (run_real_data_pilot.py), byte for byte.

    Fixed split into a 797-object development split and a 1000-object evaluated panel; drop
    features constant on the development split; standardize and fit unwhitened PCA20 on
    the development split only; the evaluated panel's PCA20 vectors are the true means.
    """
    X = load_digits().data
    dev_idx, ev_idx = np.split(np.random.default_rng(DIGITS_SPLIT_SEED).permutation(len(X)),
                               [len(X) - N_DIGITS])
    keep = X[dev_idx].std(0) > 0
    Xd, Xe = X[dev_idx][:, keep], X[ev_idx][:, keep]
    loc, scale = Xd.mean(0), Xd.std(0)
    Zd = (Xd - loc) / scale
    Vt = np.linalg.svd(Zd - Zd.mean(0), full_matrices=False)[2][:D]
    mu = ((Xe - loc) / scale - Zd.mean(0)) @ Vt.T
    return mu, median_kth_neighbor_distance(mu, K)


def synthetic_panel(family, n, dev_seed, layout_seed):
    """A predecessor family: d0 from a development draw, then a fixed evaluated layout."""
    dev = gen.draw_means(family, n, np.random.default_rng(dev_seed))
    d0 = median_kth_neighbor_distance(dev, K)
    mu = gen.draw_means(family, n, np.random.default_rng(layout_seed))
    return mu, d0


def geometry(name):
    """(mu, d0) for a named Gate 2 geometry."""
    if name == "digits":
        return digits_panel()
    if name == "digits200":
        mu, d0 = digits_panel()
        return mu[:200].copy(), d0
    if name == "surface":
        return synthetic_panel("curved_surface", 1000, SURFACE_DEV_SEED, SURFACE_LAYOUT_SEED)
    if name in gen.PRIMARY:
        fi = list(gen.PRIMARY).index(name)
        return synthetic_panel(name, 1000, SYN_DEV_BASE + fi, SYN_LAYOUT_BASE + fi)
    raise KeyError(name)
