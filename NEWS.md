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
