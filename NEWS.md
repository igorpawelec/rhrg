# rhrg 0.4.0

# rhrg 0.4.1

## Fixed

* `man/grow_crowns.Rd` still documented `max_iters = 200L` after 0.3.0 changed
  the default to `NULL`. `R CMD check` reported it as a codoc mismatch and the
  R-CMD-check workflow, which treats warnings as failures, went red on all five
  platforms. Anyone reading `?grow_crowns` was told the cap was still there.

* The cause was structural, so that is what was fixed. Every file in `man/`
  carried a hand-written header instead of roxygen2's marker, so
  `roxygenise()` refused to touch any of them and the documentation had been
  maintained by hand — free to drift from the roxygen comments it claimed to
  come from, silently, except where `codoc` happens to look. `man/` is now
  generated, and `Roxygen: list(markdown = TRUE)` is set so the markdown in
  those comments renders instead of reaching the page as `[as_pixels()]`.

  The regenerated pages are strictly richer than the ones they replace (1813
  words against 1639); nothing was lost.


* **`delineate()` now works on a real CHM file with nodata.** It did not.
  `read_chm()` returns `NA` where the file declares nodata -- terra honours the
  raster's NA flag, whatever the old docstring claimed -- `smooth_chm()` then
  propagated it and *grew* the hole, and `detect_tops()` stopped with
  "missing value where TRUE/FALSE needed", which points nowhere. The
  documented one-call entry point could not process the package's own test
  data.

  `NA` is now skipped rather than propagated: window statistics are taken over
  the cells that exist, and a window holding nothing but `NA` stays `NA`.
  `.prepare_mask()` treats `NA` as outside the mask, which is where numpy
  arrives by a different route, since `NaN > x` is `FALSE`.

  On `chm_150_2023.tif` R and pyHRG now both return **255 crowns**. Neither
  managed it before: R crashed and Python silently smoothed 3 % of the raster
  to arbitrary values, because scipy's comparison filters are undefined on
  NaN rather than NaN-aware. Matching Python would have meant reproducing that
  -- the same trap that gave this package its watershed divergence -- so both
  were defined instead. pyHRG 0.5.0 carries the other half.

  The shared cross-check gained ten scenes with nodata holes, so this is
  covered rather than assumed.

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
