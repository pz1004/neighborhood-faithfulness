"""Latent families of plan Section 7.1, with the fixed orthonormal map into D=20.

The ambient map is fitted once from a declared seed and frozen before any layout is
generated, per `synthetic_generators.all_primary_families.ambient_map`.
"""
import numpy as np

D_AMBIENT = 20
_MAP_SEED = 20260918          # frozen; the ambient map must not vary across layouts


def ambient_map(latent_dim, D=D_AMBIENT):
    """Fixed orthonormal latent -> R^D embedding. Same map for every layout."""
    rng = np.random.default_rng(_MAP_SEED)
    Q, _ = np.linalg.qr(rng.normal(size=(D, D)))
    return Q[:, :latent_dim].copy()       # D x latent_dim, orthonormal columns


def folded_curve(m, rng):
    t = rng.uniform(-1, 1, m)
    return np.column_stack([np.sin(np.pi * t), np.sin(2 * np.pi * t), 0.25 * t])


def thin_bridge(m, rng):
    u = rng.choice(3, m, p=[0.45, 0.45, 0.10])
    X = np.zeros((m, 2))
    for c, ctr in enumerate([(-2.0, 0.0), (2.0, 0.0)]):
        s = u == c
        X[s] = rng.normal(ctr, 0.25, (int(s.sum()), 2))
    s = u == 2
    X[s] = np.column_stack([rng.uniform(-2, 2, int(s.sum())),
                            rng.normal(0.0, 0.05, int(s.sum()))])
    return np.column_stack([X, np.zeros(m)])


def imbalanced_mixture(m, rng):
    u = rng.choice(3, m, p=[0.80, 0.15, 0.05])
    X = np.zeros((m, 2))
    for c, (ctr, sd) in enumerate([((-2.0, 0.0), 0.45), ((2.0, 0.0), 0.25), ((0.0, 3.0), 0.10)]):
        s = u == c
        X[s] = rng.normal(ctr, sd, (int(s.sum()), 2))
    return np.column_stack([X, np.zeros(m)])


def curved_surface(m, rng):
    s = rng.uniform(-1, 1, m)
    t = rng.uniform(-1, 1, m)
    return np.column_stack([s, t, 0.5 * (s * s + t * t)])


def gaussian_control(m, rng):
    """Mandatory full-dimensional stress control (plan Section 7.1). Already in R^D."""
    return rng.normal(size=(m, D_AMBIENT))


# ---------------------------------------------------------------------------
# Extension, September 22, 2026, for `docs/replicate_budget_guideline.md`.
#
# Trigger: that note's Section 8 records its own weak link -- the guideline is indexed
# by local intrinsic dimension, but the four primary families span only 1.03 to 2.06,
# with the next measured point the degenerate control near 14.8. Nothing measured the
# range from 3 to 10, so the table licensed no interpolation across it.
#
# Two parameterised families fill it. `paraboloid(d)` is the d-dimensional
# generalisation of `curved_surface` (which is exactly d=2), so it extends an existing
# family rather than introducing an unrelated one. `uniform_ball(d)` is structureless at
# the same intrinsic dimension, and exists to test whether the guideline's index is
# really the dimension or is partly the family's shape.
#
# Latent dimensions were chosen ONLY so the measured Levina-Bickel median lands on
# targets spanning the gap; the criterion was declared before any coverage was computed
# and involves no certificate quantity. Measured (n=1000, k=20, 3 seeds):
#
#     latent d    4       6       8      11
#     paraboloid  3.72    5.36    6.93    9.04
#     ball        3.77    5.45    7.03    9.02
# ---------------------------------------------------------------------------

def paraboloid(d):
    """Uniform on [-1,1]^d lifted by 0.5*||s||^2. d=2 reproduces `curved_surface`."""
    def f(m, rng):
        s = rng.uniform(-1, 1, (m, d))
        return np.column_stack([s, 0.5 * (s * s).sum(axis=1)])
    return f


def uniform_ball(d):
    """Uniform on the d-dimensional unit ball: intrinsic dimension d, no structure."""
    def f(m, rng):
        g = rng.normal(size=(m, d))
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        return g * rng.uniform(0, 1, (m, 1)) ** (1.0 / d)
    return f


BRIDGE_DIMS = (4, 6, 8, 11)          # frozen; see the table above
BRIDGE = {f"paraboloid_d{d}": paraboloid(d) for d in BRIDGE_DIMS}
BRIDGE.update({f"ball_d{d}": uniform_ball(d) for d in BRIDGE_DIMS})

PRIMARY = {
    "folded_curve": folded_curve,
    "thin_bridge": thin_bridge,
    "imbalanced_mixture": imbalanced_mixture,
    "curved_surface": curved_surface,
}
STRESS = {"gaussian_control": gaussian_control}
ALL = {**PRIMARY, **STRESS}
ALL_WITH_BRIDGE = {**PRIMARY, **BRIDGE, **STRESS}


def latent_dim(family):
    """Latent dimension of a family's generator, before the ambient map."""
    return ALL_WITH_BRIDGE[family](2, np.random.default_rng(0)).shape[1]


def draw_means(family, m, rng):
    """Clean feature means mu in R^20 for one layout."""
    Z = ALL_WITH_BRIDGE[family](m, rng)
    if Z.shape[1] == D_AMBIENT:
        return Z                                   # stress control is already ambient
    return Z @ ambient_map(Z.shape[1]).T
