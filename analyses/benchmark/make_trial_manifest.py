#!/usr/bin/env python3
"""Trial manifest and mouse overlap accounting (docs/revision_protocol.md, Section 6, items 1 and 3).

EEG (no data value is read: the positions, sessions and blocks follow from the published code):
- every stored repetition's merged index and session, replayed as reach/things.py does;
- its chronological rank within its session, replayed from the dataset authors' epoching step;
- its block, and its place in block U, in each of the 10 confirmatory splits.
The authors' preprocessing code (github.com/gifale95/eeg_encoding, 02_eeg_preprocessing) selects
each test image's repetitions within a session with
`idx = shuffle(np.where(events == img_cond[i])[0], random_state=seed, n_samples=20)`, seed 20200220,
the same seed as the cross-session shuffle. With exactly 20 presentations per image and session,
as the dataset documents, the same permutation applies to every image. Stored repetition j of a
session is then its chronological occurrence perm[j], with
perm = shuffle(arange(20), random_state=20200220).
Onset times are not recoverable from the preprocessed files.

Mouse V1 (the NWB files are read, each checked against its SHA-256): every presentation's image,
rank in time within its image, onset frame, epoch (the three runs of natural scenes) and stratum,
and its block in each split. Then, per experiment, the overlap of consecutive response windows,
and per split, how many consecutive presentation pairs fall one in U and one in H.

Writes results/manifest/trials_eeg.csv, trials_mouse.csv, mouse_overlap.json and overlap.txt.
"""
import csv, hashlib, json, os
import numpy as np
import h5py
from sklearn.utils import shuffle as sk_shuffle

from reach import things, allen
import run_things_confirmatory as RT
import run_allen_confirmatory as RA

OUT = "results/manifest"
N_SPLITS = 10
GAP_FRAMES = 100                      # a gap between presentations longer than this starts a new epoch


def eeg_manifest():
    merged = sk_shuffle(np.arange(80), random_state=things.AUTHOR_SEED)
    perm20 = sk_shuffle(np.arange(20), random_state=things.AUTHOR_SEED, n_samples=20)
    ses = things.sessions()
    assert (merged // 20 == ses).all()
    rows = []
    for s in range(1, 11):
        blocks = [things.split_blocks(ses, np.random.default_rng(RT.SPLIT_BASE + 1000 * s + j)) for j in range(N_SPLITS)]
        for p in range(80):
            r = dict(subject=s, position=p, merged_index=int(merged[p]), session=int(ses[p]) + 1,
                     within_session_index=int(merged[p] % 20), chronological_rank_in_session=int(perm20[merged[p] % 20]))
            for j, b in enumerate(blocks):
                for name in ("F", "U", "H"):
                    hit = np.flatnonzero(b[name] == p)
                    if len(hit):
                        r[f"split{j}"] = name
                        r[f"split{j}_order"] = int(hit[0])        # place within the block; U's first R are the analyst's
            rows.append(r)
    return rows, perm20


def presentations(path):
    """Non-blank and blank presentations of an NWB file: image (-1 blank) and start frame, in time order."""
    with h5py.File(path, "r") as f:
        img = f[allen.NS + "/data"][()].astype(int)
        start = f[allen.NS + "/frame_duration"][()][:, 0].astype(int)
        ts = f[allen.DFF + "/timestamps"][()]
    o = np.argsort(start, kind="stable")
    return img[o], start[o], 1.0 / np.median(np.diff(ts))


def verify(entry):
    path = os.path.join(RA.DATA, f"{entry['experiment_id']}.nwb")
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    assert h.hexdigest() == entry["sha256"], f"{entry['experiment_id']}: SHA-256 mismatch"
    return path


def mouse_manifest():
    rows, overlap = [], []
    st = allen.strata()
    for s, entry in enumerate(RA.included(), start=1):
        img, start, fps = presentations(verify(entry))
        epoch = np.concatenate([[1], 1 + np.cumsum(np.diff(start) > GAP_FRAMES)])
        rank = np.full(len(img), -1)
        for i in range(allen.N_IMAGES):
            sel = np.flatnonzero(img == i)                  # already in time order
            assert len(sel) == allen.N_REPS
            rank[sel] = np.arange(allen.N_REPS)
        blocks = [allen.split_blocks(st, np.random.default_rng(RA.SPLIT_BASE + 1000 * s + j)) for j in range(N_SPLITS)]
        where = [{int(p): name for name in ("F", "U", "H") for p in b[name]} for b in blocks]
        for t in np.flatnonzero(img >= 0):
            r = dict(subject=s, experiment_id=entry["experiment_id"], image=int(img[t]), rank_in_image=int(rank[t]),
                     onset_frame=int(start[t]), epoch=int(epoch[t]), stratum="early" if st[rank[t]] == 0 else "late")
            for j in range(N_SPLITS):
                r[f"split{j}"] = where[j].get(int(rank[t]), "unused")
            rows.append(r)
        # consecutive presentations: window overlap, and pairs that fall one in U and one in H
        nwin = {w: int(round(sec * fps)) for w, sec in allen.WINDOWS.items()}
        step = np.diff(start)
        both = (img[:-1] >= 0) & (img[1:] >= 0)
        o = dict(subject=s, experiment_id=entry["experiment_id"], fps=float(fps), n_presentations=int(len(img)),
                 blanks=int((img < 0).sum()), epochs=int(epoch.max()), median_step_frames=float(np.median(step)),
                 window_frames=nwin,
                 overlap_frames={w: float(np.median(np.maximum(0, n - step[both]))) for w, n in nwin.items()},
                 consecutive_pairs=int(both.sum()), same_image_pairs=int((both & (img[:-1] == img[1:])).sum()))
        uh, uh_same = [], []
        for j in range(N_SPLITS):
            blk = np.array([where[j].get(int(rk), "unused") if im >= 0 else "blank" for im, rk in zip(img, rank)])
            cross = both & (((blk[:-1] == "U") & (blk[1:] == "H")) | ((blk[:-1] == "H") & (blk[1:] == "U")))
            uh.append(int(cross.sum()))
            uh_same.append(int((cross & (img[:-1] == img[1:])).sum()))
        n_pairs = allen.N_IMAGES * (allen.N_IMAGES - 1) / 2
        o.update(u_h_pairs=uh, u_h_share=float(np.mean(uh) / both.sum()), u_h_same_image=uh_same,
                 u_h_per_image_pair=float(np.mean(uh) / n_pairs))
        overlap.append(o)
    return rows, overlap


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    eeg, perm20 = eeg_manifest()
    write_csv(f"{OUT}/trials_eeg.csv", eeg)
    mouse, overlap = mouse_manifest()
    write_csv(f"{OUT}/trials_mouse.csv", mouse)
    json.dump(overlap, open(f"{OUT}/mouse_overlap.json", "w"), indent=1)
    lines = ["TRIAL MANIFEST AND OVERLAP ACCOUNTING (make_trial_manifest.py)",
             f"EEG: {len(eeg)} rows (10 subjects x 80 stored repetitions); within-session permutation of the"
             f" authors' epoching step: {perm20.tolist()}",
             f"Mouse V1: {len(mouse)} rows (10 experiments x 118 images x 50 presentations)",
             "Mouse V1, per experiment: frames per window (0.5 s / 0.25 s), median overlap of consecutive windows,"
             " consecutive pairs one in U and one in H (mean over the 10 splits), and per image pair"]
    for o in overlap:
        lines.append(f"  subject {o['subject']:2d}: fps {o['fps']:.2f}, step {o['median_step_frames']:.0f} frames,"
                     f" windows {o['window_frames']['full']}/{o['window_frames']['early']} frames,"
                     f" overlap {o['overlap_frames']['full']:.0f}/{o['overlap_frames']['early']:.0f};"
                     f" {o['epochs']} epochs; U-H pairs {np.mean(o['u_h_pairs']):.0f} of {o['consecutive_pairs']}"
                     f" ({o['u_h_share']:.3f}), {o['u_h_per_image_pair']:.3f} per image pair,"
                     f" same image {np.mean(o['u_h_same_image']):.1f}")
    open(f"{OUT}/overlap.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
