"""THINGS-EEG2 preprocessed test set for the need pilot (plan Amendment 2).

Safe loading, reconstruction of session membership, session-stratified blocks, and the
representation fitted on block F only.
"""
import pickle
import numpy as np
from numpy.lib import format as npf
from sklearn.utils import shuffle as sk_shuffle
from sklearn.covariance import ledoit_wolf

AUTHOR_SEED = 20200220                 # seed of the authors' cross-session shuffle
N_SES, REPS_PER_SES = 4, 20
BLOCKS = (("F", 2), ("U", 8), ("H", 10))   # repetitions per session in each block

# The file is np.save of a dict, i.e. a pickle. Only NumPy's array reconstruction is let
# through; any other global in the stream raises instead of executing.
_ALLOWED = {("numpy.core.multiarray", "_reconstruct"), ("numpy._core.multiarray", "_reconstruct"),
            ("numpy.core.multiarray", "scalar"), ("numpy._core.multiarray", "scalar"),
            ("numpy", "ndarray"), ("numpy", "dtype")}


class _Restricted(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) not in _ALLOWED:
            raise pickle.UnpicklingError(f"blocked global {module}.{name}")
        if name in ("_reconstruct", "scalar"):
            return getattr(np._core.multiarray, name)
        return getattr(np, name)


def load_test(path):
    """(X, ch_names, times) with X of shape images x repetitions x channels x samples."""
    with open(path, "rb") as f:
        version = npf.read_magic(f)
        read = npf.read_array_header_1_0 if version == (1, 0) else npf.read_array_header_2_0
        _, _, dtype = read(f)
        if not dtype.hasobject:
            raise ValueError("expected a pickled object array")
        obj = _Restricted(f).load()
    d = obj.item() if isinstance(obj, np.ndarray) else obj
    return np.asarray(d["preprocessed_eeg_data"]), list(d["ch_names"]), np.asarray(d["times"])


def sessions():
    """Session of every stored repetition, by replaying the authors' shuffle.

    Their save step merges sessions along the repetition axis (session 1 first) and then
    applies merged[:, idx] with idx = sklearn.utils.shuffle(arange(80), random_state=seed),
    so stored position p holds merged repetition idx[p].
    """
    idx = sk_shuffle(np.arange(N_SES * REPS_PER_SES), random_state=AUTHOR_SEED)
    return idx // REPS_PER_SES


def split_blocks(ses, rng):
    """Positions of blocks F, U, H, each ordered so that every prefix of 4j is stratified.

    Within each session the repetitions are permuted and dealt to F, U and H. A block's
    order is round-robin over sessions, so its first 4j entries hold j per session.
    """
    per = {name: [] for name, _ in BLOCKS}
    for s in range(N_SES):
        pos = rng.permutation(np.flatnonzero(ses == s))
        o = 0
        for name, m in BLOCKS:
            per[name].append(pos[o:o + m])
            o += m
    return {name: np.stack(per[name], axis=1).reshape(-1) for name, _ in BLOCKS}


def features(X, times):
    """Post-stimulus samples (t >= 0), channels x samples flattened: images x reps x p."""
    t0 = int(np.flatnonzero(times >= -1e-9)[0])
    return X[:, :, :, t0:].reshape(X.shape[0], X.shape[1], -1)


def fit_representation(XF, d_pca=100, D=20):
    """Linear map to D dimensions fitted on block F only (plan Amendment 2, item 2).

    PCA of F's image means to d_pca dimensions; F's residual covariance there (Ledoit-Wolf,
    rescaled to its n(r-1) degrees of freedom); whitening; PCA of the whitened means to D.
    """
    n, r, p = XF.shape
    M = XF.mean(axis=1)
    m0 = M.mean(axis=0)
    P = np.linalg.svd(M - m0, full_matrices=False)[2][:d_pca].T
    E = (XF - M[:, None, :]).reshape(-1, p) @ P
    cov, shrink = ledoit_wolf(E, assume_centered=True)
    cov *= (n * r) / (n * (r - 1))
    lam, vec = np.linalg.eigh(cov)
    Wh = (vec / np.sqrt(lam)) @ vec.T
    Z = (M - m0) @ P @ Wh
    V = np.linalg.svd(Z - Z.mean(axis=0), full_matrices=False)[2][:D].T
    return P @ Wh @ V, m0, float(shrink)
