#!/usr/bin/env python3
"""Describe the 10 included Allen experiments for the manuscript's Methods: donor, sex, depth and
cell count, from the downloaded experiment metadata (not tracked) and the tracked manifest.

Writes results/allen_confirmatory/subjects_described.json, in the analysis's subject order.
"""
import json

subs = json.load(open("results/allen_confirmatory/subjects.json"))
man = {m["experiment_id"]: m for m in json.load(open("data/allen_2p/manifest.json"))}
meta = {int(e["id"]): e for e in json.load(open("data/allen_2p/ophys_experiments.json"))}


def donor_sex(e):
    """The first 'sex' field found in the experiment's nested specimen record."""
    stack = [e]
    while stack:
        d = stack.pop()
        if isinstance(d, dict):
            if isinstance(d.get("sex"), str):
                return d["sex"]
            stack += list(d.values())
        elif isinstance(d, list):
            stack += d
    raise KeyError("no sex field")


out = []
for s in subs:
    e, m = meta[int(s["experiment_id"])], man[s["experiment_id"]]
    out.append(dict(subject=s["subject"], experiment_id=s["experiment_id"], donor_id=m["donor_id"],
                    sex=donor_sex(e), imaging_depth=e["imaging_depth"], n_cells=s["n_cells"]))
json.dump(out, open("results/allen_confirmatory/subjects_described.json", "w"), indent=1)
for r in out:
    print(r)
