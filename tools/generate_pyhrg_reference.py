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
# Cases 60-69 are large and crowded on purpose. The original 60 are 30-60 px
# with 2-7 trees, which is nothing like a real canopy height model: the
# divergence found on real CHMs in 2026-07 lives in plateau structure that
# only appears once crowns are dense enough to touch.
for k in range(70):
    big = k >= 60
    n = int(rng.integers(140, 200)) if big else int(rng.integers(30, 60))
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))
    n_trees = int(rng.integers(25, 60)) if big else int(rng.integers(2, 8))
    for _ in range(n_trees):
        r, c = rng.integers(4, n - 4, 2)
        chm = np.maximum(chm, rng.uniform(10, 28) *
                         np.exp(-((yy - r) ** 2 + (xx - c) ** 2)
                                / (2 * rng.uniform(3, 7) ** 2)))
    if k % 3 == 0:
        chm = np.round(chm * 2) / 2          # plateaus, where ties are common
    if k % 7 == 0:
        # Nodata holes. Until pyHRG 0.5.0 / rHRG 0.4.0 these two disagreed
        # completely here: scipy's comparison filters are undefined on NaN,
        # so Python smoothed the hole edges to arbitrary window elements,
        # while R propagated NA and then stopped with "missing value where
        # TRUE/FALSE needed". Both now skip NaN, and this is what checks it.
        chm[2:6, 2:6] = np.nan
        chm[n - 7:n - 3, n - 9:n - 5] = np.nan
    chm = chm.astype(np.float32)

    sm = smooth_chm(chm, ws=3, method="median")
    sub = detect_tops(sm, hmin=5, ws=5)
    tops = as_pixels(sub)
    if not tops:
        continue
    vt = float(rng.choice([2.0, 8.0, 25.0]))
    rule = str(rng.choice(["height", "distance", "similarity"]))

    g = HierarchicalRegionGrower(sm)
    out = g.run_all(tops, variance_thresh=vt, mask_thresh=1.0,
                    conflict_rule=rule)
    np.savetxt(f"{OUT}/sm_{k}.csv", sm, delimiter=",")
    np.savetxt(f"{OUT}/tops_{k}.csv", np.array(tops) + 1, fmt="%d", delimiter=",")
    # Subpixel tops, before as_pixels. Without these the R side of
    # detect_tops is never compared against anything: the validator used to
    # read tops that Python had already detected *and* floored, so R's own
    # detect_tops -> grow_crowns handoff went unexercised.
    np.savetxt(f"{OUT}/sub_{k}.csv", np.asarray(sub, float) + 1.0,
               delimiter=",", fmt="%.17g")
    np.savetxt(f"{OUT}/out_{k}.csv", out, fmt="%d", delimiter=",")
    cases.append((k, vt, rule, g.n_contested))

with open(f"{OUT}/meta.csv", "w") as f:
    for k, vt, rule, nc in cases:
        f.write(f"{k},{vt},{rule},{nc}\n")
print(f"{len(cases)} reference cases written to {OUT}/")
