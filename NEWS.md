# rhrg 0.3.0

* **`max_iters` now defaults to `NULL`, meaning natural termination.** It was
  200, and the bound bit: on `chm_33_2012.tif` 332 of 492 crowns stopped there
  with candidates still queued, and the crown count read **132 against 63**
  once the cap was lifted -- more than a factor of two in the headline number,
  decided by a constant rather than by the canopy.

  The boundaries barely moved, 2.9 % of the partition, because a truncated
  grow blocks merges rather than misplacing pixels. That is what made it hard
  to notice: the segmentation looked right and the tree count did not.

  The cap bought nothing. Growth is bounded anyway -- each iteration either
  accepts a region or records a rejection, and there are finitely many of both
  -- and natural termination needed at most 484 iterations on that scene while
  running *faster*, since twice as many surviving crowns cost more in conflict
  arbitration than the extra merges cost in growing.

  Passing an explicit `max_iters` still works and now warns when it binds.

  pyHRG 0.3.0 carries the same change.

# rhrg 0.2.0

* **An even smoothing or detection window is now refused.** An even window has
  no centre pixel, so it sits half a pixel off and the result depends on which
  way the raster happens to be oriented. Measured before the guard: smoothing a
  40x55 scene and its mirror image differed by up to 10.5 m at `ws = 4`, and
  `detect_tops()` found 151 tops against 137 on the mirror. `method =
  "gaussian"` is exempt, because `ws` only sets sigma there and the kernel stays
  symmetric.

  pyHRG 0.2.0 carries the same guard; the asymmetry was faithfully reproduced
  here, so it had to be fixed on both sides to keep the two in step.

# rhrg 0.1.0

First release.

* `delineate()` runs the whole pipeline; `smooth_chm()`, `detect_tops()`,
  `merge_tops()`, `screen_tops()` and `grow_crowns()` are each usable alone.
* `conflict_rule` decides which crown wins canopy claimed by two: `"height"`,
  `"distance"` or `"similarity"`. Ties resolve to the lower crown id, so the
  segmentation does not depend on the order the tree tops arrive in.
* `protect_seeds` guarantees one crown per tree top, disabling merging.
* `retry_rejected` reconsiders regions rejected earlier in a grow.
* `read_chm()`, `write_crowns()` and `crowns_to_polygons()` bridge to
  \pkg{terra}, which is only a Suggests: the algorithm works on plain matrices
  with no spatial dependency.
* Verified identical to the pyHRG Python package on 60 random scenes
  (120,977 pixels, no difference), and the watershed alone on a further 200
  scenes against `scikit-image`. See `tools/cross_validate_against_pyhrg.R`.
* On real canopy height models the agreement is close but not exact: 375 of
  147,368 masked pixels differ in the watershed (0.25 %), which the region
  growing amplifies to 5,142 of 281,602 crown pixels. `src/watershed.c`
  reimplements `skimage.segmentation.watershed` rather than calling it, and
  the two break plateau ties differently. Documented in the README rather
  than described as equality.
* `tools/measure_tie_sensitivity.py` quantifies how much the method itself
  depends on those ties. A nanometre of noise removes 3 to 30 tree tops and
  moves crowns for up to 62 % of canopy pixels — but 0 % once the tops are
  held fixed, so the instability is in tree detection, not in the growing.
