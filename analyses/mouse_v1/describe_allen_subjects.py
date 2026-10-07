#!/usr/bin/env python3
"""Describe the 10 included Allen experiments for the manuscript's Methods: donor, sex, genotype,
age at the recording session, depth and cell count, from the downloaded experiment metadata (not
tracked) and the tracked manifest. Sex and genotype are the donor record's fields; the age is the
number of days from the donor's date of birth to the experiment's date of acquisition.

Writes results/allen_confirmatory/subjects_described.json, in the analysis's subject order.
"""
import datetime as dt
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


def day(stamp):
    return dt.datetime.fromisoformat(stamp.replace("Z", "+00:00"))


out = []
for s in subs:
    e, m = meta[int(s["experiment_id"])], man[s["experiment_id"]]
    donor = e["specimen"]["donor"]
    assert donor["id"] == m["donor_id"] and donor["sex"] == donor_sex(e)
    out.append(dict(subject=s["subject"], experiment_id=s["experiment_id"], donor_id=m["donor_id"],
                    sex=donor_sex(e), genotype=donor["full_genotype"],
                    age_days=(day(e["date_of_acquisition"]) - day(donor["date_of_birth"])).days,
                    imaging_depth=e["imaging_depth"], n_cells=s["n_cells"]))
json.dump(out, open("results/allen_confirmatory/subjects_described.json", "w"), indent=1)
for r in out:
    print(r)
