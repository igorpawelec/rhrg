"""Write reference crowns from pyHRG for tools/cross_validate_against_pyhrg.R.

    pip install pyhrg
    python3 tools/generate_pyhrg_reference.py

Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
"""
import os
import numpy as np
from pyhrg import smooth_chm, detect_tops, as_pixels, HierarchicalRegionGrower

OUT = os.path.join("tools", "reference")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(11)
cases = []
for k in range(60):
    n = int(rng.integers(30, 60))
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))
    for _ in range(int(rng.integers(2, 8))):
        r, c = rng.integers(4, n - 4, 2)
        chm = np.maximum(chm, rng.uniform(10, 28) *
                         np.exp(-((yy - r) ** 2 + (xx - c) ** 2)
                                / (2 * rng.uniform(3, 7) ** 2)))
    if k % 3 == 0:
        chm = np.round(chm * 2) / 2          # plateaus, where ties are common
    chm = chm.astype(np.float32)

    sm = smooth_chm(chm, ws=3, method="median")
    tops = as_pixels(detect_tops(sm, hmin=5, ws=5))
    if not tops:
        continue
    vt = float(rng.choice([2.0, 8.0, 25.0]))
    rule = str(rng.choice(["height", "distance", "similarity"]))

    g = HierarchicalRegionGrower(sm)
    out = g.run_all(tops, variance_thresh=vt, mask_thresh=1.0,
                    conflict_rule=rule)
    np.savetxt(f"{OUT}/sm_{k}.csv", sm, delimiter=",")
    np.savetxt(f"{OUT}/tops_{k}.csv", np.array(tops) + 1, fmt="%d", delimiter=",")
    np.savetxt(f"{OUT}/out_{k}.csv", out, fmt="%d", delimiter=",")
    cases.append((k, vt, rule, g.n_contested))

with open(f"{OUT}/meta.csv", "w") as f:
    for k, vt, rule, nc in cases:
        f.write(f"{k},{vt},{rule},{nc}\n")
print(f"{len(cases)} reference cases written to {OUT}/")
