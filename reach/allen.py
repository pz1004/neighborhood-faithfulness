"""Allen Brain Observatory two-photon natural scenes (docs/note_protocol_allen.md).

Per-presentation responses from an NWB 1.0 file, the two time strata, stratum-balanced
blocks and cross-validation halves.
"""
import numpy as np
import h5py

N_IMAGES, N_REPS = 118, 50
N_STRATA = 2                                   # early and late half of each image's presentations
BLOCKS = (("F", 5), ("U", 8), ("H", 12))       # repetitions per stratum in each block
WINDOWS = {"full": 0.5, "early": 0.25}         # response windows after onset, seconds

DFF = "processing/brain_observatory_pipeline/DfOverF/imaging_plane_1"
NS = "stimulus/presentation/natural_scenes_stimulus"


def load_responses(path, windows=WINDOWS):
    """{window: X} with X of shape images x repetitions x cells, repetitions in time order.

    A presentation's response is each cell's mean dF/F over the imaging frames from its
    first frame (frame_duration[:, 0]) for the window's length at the recording's frame rate.
    Blank presentations (image index -1) are dropped.
    """
    with h5py.File(path, "r") as f:
        dff = f[DFF + "/data"][()]                      # cells x frames
        ts = f[DFF + "/timestamps"][()]
        img = f[NS + "/data"][()].astype(int)
        start = f[NS + "/frame_duration"][()][:, 0].astype(int)
    fps = 1.0 / np.median(np.diff(ts))
    keep = img >= 0
    img, start = img[keep], start[keep]
    out = {}
    for name, sec in windows.items():
        nwin = int(round(sec * fps))
        X = np.empty((N_IMAGES, N_REPS, dff.shape[0]))
        for i in range(N_IMAGES):
            onsets = np.sort(start[img == i])           # time order within image
            if len(onsets) != N_REPS:
                raise ValueError(f"image {i}: {len(onsets)} presentations, expected {N_REPS}")
            for r, s in enumerate(onsets):
                X[i, r] = dff[:, s:s + nwin].mean(axis=1)
        out[name] = X
    return out, dict(fps=float(fps), n_cells=int(dff.shape[0]))


def strata():
    """Stratum of each repetition position: the first 25 in time, then the last 25."""
    return np.repeat(np.arange(N_STRATA), N_REPS // N_STRATA)


def split_blocks(st, rng):
    """Positions of blocks F, U, H; each ordered round-robin over the strata."""
    per = {name: [] for name, _ in BLOCKS}
    for s in range(N_STRATA):
        pos = rng.permutation(np.flatnonzero(st == s))
        o = 0
        for name, m in BLOCKS:
            per[name].append(pos[o:o + m])
            o += m
    return {name: np.stack(per[name], axis=1).reshape(-1) for name, _ in BLOCKS}


def halves(R):
    """Two stratum-balanced halves of a round-robin prefix of length R (R % 4 == 0)."""
    p = np.arange(R)
    return p[(p // N_STRATA) % 2 == 0], p[(p // N_STRATA) % 2 == 1]
