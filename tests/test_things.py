#!/usr/bin/env python3
"""Checks of reach/things.py that touch no THINGS-EEG2 data value (plan Amendment 2)."""
import io, os, pickle, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from reach import things

fails = []


def check(name, cond):
    print(f"[{'ok' if cond else 'FAIL'}] {name}")
    if not cond:
        fails.append(name)


with tempfile.TemporaryDirectory() as d:
    # a file with the authors' layout, written the way np.save writes a dict
    rng = np.random.default_rng(3_090_001)
    X = rng.normal(size=(5, 80, 17, 100))
    p = os.path.join(d, "ok.npy")
    np.save(p, {"preprocessed_eeg_data": X, "ch_names": ["O1", "Oz"], "times": np.arange(-20, 80) / 100})
    Y, ch, t = things.load_test(p)
    check("restricted loader reads an np.save'd dict exactly", np.array_equal(Y, X) and ch == ["O1", "Oz"])

    # a hostile file: an object array whose pickle calls os.system
    class Evil:
        def __reduce__(self):
            return (os.system, ("echo pwned",))
    q = os.path.join(d, "evil.npy")
    np.save(q, np.array([Evil()], dtype=object), allow_pickle=True)
    try:
        things.load_test(q)
        check("restricted loader blocks non-numpy globals", False)
    except pickle.UnpicklingError:
        check("restricted loader blocks non-numpy globals", True)

ses = things.sessions()
check("replayed shuffle gives 4 sessions of 20", sorted(np.bincount(ses).tolist()) == [20, 20, 20, 20])
b = things.split_blocks(ses, np.random.default_rng(0))
allpos = np.concatenate([b["F"], b["U"], b["H"]])
check("blocks partition the 80 repetitions", sorted(allpos.tolist()) == list(range(80)))
strat = all((np.bincount(ses[b["U"][:4 * j]], minlength=4) == j).all() for j in (2, 4, 8))
strat &= all((np.bincount(ses[b["H"][j:j + Rh]], minlength=4) == Rh // 4).all()
             for Rh in (8, 20, 40) for j in range(0, 40, Rh))
check("user prefixes and reference sub-blocks are session-stratified", strat)
t = np.arange(-20, 80) / 100
check("features keep t >= 0 only (80 samples x 17 channels)",
      things.features(np.zeros((2, 3, 17, 100)), t).shape == (2, 3, 1360))
XF = rng.normal(size=(200, 8, 60)) + rng.normal(size=(200, 1, 60)) * 3
W, m0, sh = things.fit_representation(XF, d_pca=30, D=5)
E = ((XF - XF.mean(1, keepdims=True)) @ W).reshape(-1, 5)
check("representation whitens F's own noise (variance ~ 1)", abs(E.var(0).mean() * 8 / 7 - 1) < 0.1)
print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
