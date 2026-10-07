#!/usr/bin/env python3
"""Select and fetch the Allen two-photon recordings of docs/note_protocol_allen.md, Section 2.

Candidates come from the experiment metadata (data/allen_2p/ophys_experiments.json):
- session B, VISp, Cux2-CreERT2, imaging depth 175 or 275 um;
- one experiment per mouse, the 175 um one if the mouse has both;
- ordered by ascending donor ID.

Files are downloaded in that order from Allen's public bucket. Each one passes the
inclusion check, which reads only the ROI count and the natural-scenes stimulus table,
never a response value. Downloading stops when 10 experiments are included. SHA-256s and
inclusion results go to data/allen_2p/manifest.json.

    python3 analyses/fetch/fetch_allen_2p.py --list        # candidates and file sizes only
    python3 analyses/fetch/fetch_allen_2p.py               # download until 10 are included
"""
import argparse, hashlib, json, os, re, sys, time, urllib.request

BUCKET = "https://allen-brain-observatory.s3.us-west-2.amazonaws.com"
META = "data/allen_2p/ophys_experiments.json"
OUT = "data/allen_2p"
N_INCLUDE, MIN_CELLS, N_IMAGES, N_REPS = 10, 100, 118, 50
BYTE_CAP = 5.0e9            # stop and ask the owner if inclusion needs more than this


def cre(x):
    ls = [t["name"] for t in x["specimen"]["donor"].get("transgenic_lines", [])
          if t.get("transgenic_line_type_name") == "driver" and "Cre" in t["name"]]
    return ls[0] if ls else None


def candidates():
    d = json.load(open(META))
    pool = [x for x in d if x["stimulus_name"] == "three_session_B"
            and x["targeted_structure"]["acronym"] == "VISp"
            and cre(x) == "Cux2-CreERT2" and x["imaging_depth"] in (175, 275)]
    per_mouse = {}
    for x in sorted(pool, key=lambda x: (x["imaging_depth"], x["id"])):   # 175 um first
        per_mouse.setdefault(x["specimen"]["donor_id"], x)
    return [per_mouse[k] for k in sorted(per_mouse)]


def s3_size(key):
    req = urllib.request.Request(f"{BUCKET}/{key}", method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def inclusion(path):
    """(n_cells, presentations-per-image check) from the file's metadata, no responses read.

    NWB 1.0 layout of these files: natural_scenes_stimulus/data holds one image index per
    presentation (-1 = blank); frame_duration holds its first and last imaging frames.
    """
    import h5py, numpy as np
    with h5py.File(path, "r") as f:
        n_cells = int(f["processing/brain_observatory_pipeline/DfOverF/imaging_plane_1/data"].shape[0])
        img = f["stimulus/presentation/natural_scenes_stimulus/data"][()].astype(int)
    counts = np.bincount(img[img >= 0], minlength=N_IMAGES)
    ok_table = len(counts) == N_IMAGES and bool((counts == N_REPS).all())
    return n_cells, ok_table


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(key, dest):
    tmp = dest + ".part"
    h = hashlib.sha256()
    with urllib.request.urlopen(f"{BUCKET}/{key}", timeout=600) as r, open(tmp, "wb") as f:
        while True:
            chunk = r.read(8 * 1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
            f.write(chunk)
    os.replace(tmp, dest)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    cands = candidates()
    rows, total = [], 0
    for x in cands:
        key = f"visual-coding-2p/ophys_experiment_data/{x['id']}.nwb"
        size = s3_size(key)
        total += size
        rows.append(dict(experiment_id=x["id"], donor_id=x["specimen"]["donor_id"],
                         depth=x["imaging_depth"], container_id=x["experiment_container_id"],
                         key=key, bytes=size))
        print(f"donor {x['specimen']['donor_id']}  experiment {x['id']}  depth {x['imaging_depth']} um  "
              f"{size / 1e6:6.1f} MB", flush=True)
    print(f"{len(rows)} candidate mice; all files together {total / 1e9:.2f} GB")
    if a.list:
        return
    manifest, included, spent = [], 0, 0
    for r in rows:
        if included == N_INCLUDE:
            break
        if spent + r["bytes"] > BYTE_CAP:
            print(f"byte cap reached ({spent / 1e9:.2f} GB so far): stopping to ask the owner")
            break
        spent += r["bytes"]
        dest = os.path.join(OUT, f"{r['experiment_id']}.nwb")
        t0 = time.time()
        sha = sha256_file(dest) if os.path.exists(dest) else fetch(r["key"], dest)
        n_cells, ok_table = inclusion(dest)
        inc = n_cells >= MIN_CELLS and ok_table
        included += inc
        manifest.append(dict(r, sha256=sha, n_cells=n_cells, table_50x118=ok_table, included=inc,
                             fetched=time.strftime("%Y-%m-%d %H:%M:%S")))
        print(f"  {r['experiment_id']}: {r['bytes'] / 1e6:.0f} MB in {time.time() - t0:.0f}s, "
              f"{n_cells} cells, 50x118 table {ok_table} -> {'included' if inc else 'excluded'}", flush=True)
        json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print(f"{included} included of {len(manifest)} downloaded")


if __name__ == "__main__":
    main()
