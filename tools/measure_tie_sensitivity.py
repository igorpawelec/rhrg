"""How much do the crowns move when only the plateau ties change?

A canopy height model is full of plateaus — flat roofs of smoothed crowns,
and whole regions of ground clipped to the same value. Marker-based
watershed has to break those ties somehow, and the choice is arbitrary:
scikit-image floods first-in-first-out by insertion age, a different
implementation may not, and neither is more correct than the other.

This script puts a number on what that arbitrariness costs. It adds noise
of 1e-9 m — a nanometre, far below any physical or instrumental meaning,
and far below the float32 the data is stored in — which changes nothing
about the canopy but does break the ties differently each time. Whatever
the crowns do in response is caused by tie-breaking alone.

    python3 tools/measure_tie_sensitivity.py [n_repeats]

Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
"""
import os
import sys

import numpy as np
import rasterio
from pyhrg import (smooth_chm, detect_tops, as_pixels,
                   HierarchicalRegionGrower)

SRC = os.environ.get("CHM_DIR", r"D:\Apps\pyhrg\test_data")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
EPS = 1e-9

HMIN, WS, SMOOTH_WS = 5.0, 5, 3
VT, MASK, RULE = 8.0, 1.0, "height"


def run(sm):
    tops = as_pixels(detect_tops(sm, hmin=HMIN, ws=WS))
    if not tops:
        return None, None, 0
    g = HierarchicalRegionGrower(sm)
    out = g.run_all(tops, variance_thresh=VT, mask_thresh=MASK,
                    conflict_rule=RULE)
    return out, tops, len(tops)


def seed_map(out, tops):
    """Crown ids -> the tree top each pixel belongs to.

    Crown ids are arbitrary: id 7 in one run and id 7 in another need not
    be the same tree. Comparing the integers directly counts a renumbering
    as a total change and reports ~100% wherever the top ordering shifts,
    which is a measurement artefact, not a result. Resolving each pixel to
    its seed *coordinates* compares the partition itself.
    """
    coords = np.zeros((len(tops) + 1, 2), dtype=np.int64)
    coords[0] = (-1, -1)                       # background
    for i, (r, c) in enumerate(tops, start=1):
        coords[i] = (r, c)
    idx = np.clip(out, 0, len(tops))
    return coords[idx]


print("noise = %g m, %d repeats per scene\n" % (EPS, N))
print("%-16s %6s %10s %9s %17s %17s" %
      ("case", "tops", "canopy px", "tops mvd",
       "tops free: worst", "tops fixed: worst"))
print("-" * 84)

rows = []
for f in sorted(os.listdir(SRC)):
    if not f.endswith(".tif"):
        continue
    name = f[:-4]
    with rasterio.open(os.path.join(SRC, f)) as src:
        raw = src.read(1).astype(np.float64)
        nod = src.nodata
    if nod is not None:
        raw[raw == nod] = 0.0
    raw[~np.isfinite(raw)] = 0.0

    base_sm = smooth_chm(raw, ws=SMOOTH_WS, method="median")
    base, base_tops, n_tops = run(base_sm)
    if base is None:
        continue
    canopy = int((base > 0).sum())
    base_seed = seed_map(base, base_tops)

    diffs, fixed_diffs, top_moves = [], [], []
    for k in range(N):
        rng = np.random.default_rng(1000 + k)
        sm = smooth_chm(raw + rng.uniform(0, EPS, raw.shape),
                        ws=SMOOTH_WS, method="median")

        # (a) whole pipeline: the tops are free to move too, because a local
        #     maximum on a plateau is itself a tie.
        out, tops, _ = run(sm)
        if out is None:
            continue
        top_moves.append(len(set(tops) ^ set(base_tops)))
        moved = (seed_map(out, tops) != base_seed).any(axis=-1)
        # Only canopy counts: background staying background is not news.
        diffs.append(int((moved & ((out > 0) | (base > 0))).sum()))

        # (b) tops frozen at the baseline, so detection ties are removed and
        #     what is left is the watershed and the region growing alone.
        #     Without this split the two effects are reported as one and the
        #     watershed gets blamed for what tree detection did.
        g = HierarchicalRegionGrower(sm)
        out_f = g.run_all(base_tops, variance_thresh=VT, mask_thresh=MASK,
                          conflict_rule=RULE)
        moved_f = (seed_map(out_f, base_tops) != base_seed).any(axis=-1)
        fixed_diffs.append(int((moved_f & ((out_f > 0) | (base > 0))).sum()))

    if not diffs:
        continue
    d = np.array(diffs, float)
    fd = np.array(fixed_diffs, float) if fixed_diffs else np.zeros(1)
    print("%-16s %6d %10d %9.1f %9d %6.2f%% %9d %6.2f%%" %
          (name, n_tops, canopy, float(np.mean(top_moves)),
           int(d.max()), 100 * d.max() / max(canopy, 1),
           int(fd.max()), 100 * fd.max() / max(canopy, 1)))
    rows.append((name, n_tops, canopy, d.max(), fd.max()))

if rows:
    tot_canopy = sum(r[2] for r in rows)
    tot_free = sum(r[3] for r in rows)
    tot_fixed = sum(r[4] for r in rows)
    print("-" * 84)
    print("across %d scenes, worst case of %s canopy pixels:" %
          (len(rows), format(tot_canopy, ",")))
    print("  tops free  : %d change crown (%.2f%%)" %
          (int(tot_free), 100 * tot_free / tot_canopy))
    print("  tops fixed : %d change crown (%.2f%%)" %
          (int(tot_fixed), 100 * tot_fixed / tot_canopy))
    print("\nNothing about the canopy changed. Only which pixel won a tie.")
