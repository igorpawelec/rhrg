# R HRG

<img src="https://raw.githubusercontent.com/igorpawelec/pyhrg/main/www/logoR.png" alt="R HRG logo" align="right" width="200"/>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%203.6-blue.svg)](https://www.r-project.org)

**Individual tree crown delineation from canopy height models, by Hierarchical Region Growing.**

R and plain C. No Rcpp, no compiled dependencies, and `terra` only if you want to read files.

> **Python users:** the same algorithm is in [pyHRG](https://github.com/igorpawelec/pyhrg). The two are separate repositories because their tooling and idioms do not mix — but they are not merely similar. They produce **identical output, pixel for pixel**, and that is checked rather than claimed (see [Agreement with pyHRG](#agreement-with-pyhrg)).

## Background

Delineating individual crowns from a canopy height model runs into one persistent problem: **tree tops are over-detected**. A broad crown has a ragged upper surface, so local-maximum detection finds several peaks on it. Lower the sensitivity and you start losing real trees instead.

rHRG treats that as the central problem rather than a preprocessing nuisance. Surplus tops are allowed, and the growing merges them back:

1. **Watershed.** Tree tops seed a marker-based watershed on the inverted CHM, yielding *exactly one region per tree top* — so the regions **are** the detected trees.
2. **Region adjacency graph.** Neighbouring regions get a weighted edge,
   *w(a,b) = α·|Δμ| + β·|Δσ| + γ/(border+1)*. Region statistics come from the same single pixel pass.
3. **Growing.** Each tree grows greedily, absorbing neighbours while the combined height variance stays under `variance_thresh`. When two trees absorb each other, they are one crown — this is how over-detection is corrected.
4. **Arbitration.** Growing runs independently per seed, so two crowns can claim the same region. Those claims are settled on the data (taller tree, nearest tree, or best height match), not on iteration order, so the result does not depend on how the tree tops were sorted.

`variance_thresh` is the main control: it sets how much height variation one crown may contain, and therefore how readily neighbouring trees merge.

## Installation

```r
# install.packages("remotes")
remotes::install_github("igorpawelec/rhrg")
```

A C compiler is needed, which R already requires on Linux and macOS; on Windows install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

`terra` is optional and only used to read and write rasters:

```r
install.packages("terra")
```

## Quick start

```r
library(rhrg)

res <- delineate("chm.tif",
                 smooth_ws      = 3,
                 hmin           = 7,
                 detect_ws      = 5,
                 merge_distance = 5,
                 screen_hmin    = 10,
                 variance_thresh = 2)

res$crowns        # SpatRaster of crown ids, 0 = background
res$tops          # tree tops actually used
res$n_contested   # regions claimed by more than one crown

write_crowns(res$crowns, "crowns.tif")
plot(crowns_to_polygons(res$crowns))
```

### Matrices, without touching the disk

Every stage is a plain function, and none of the algorithm needs `terra`:

```r
sm   <- smooth_chm(chm, ws = 3, method = "median")
tops <- detect_tops(sm, hmin = 7, ws = 5)      # (n, 2) subpixel row/col
tops <- merge_tops(tops, distance = 5)

crowns <- grow_crowns(sm, as_pixels(tops), variance_thresh = 2)
attr(crowns, "n_contested")
```

## Parameters

**Smoothing** — `smooth_chm(chm, ws, method)`

| Parameter | Default | Description |
|---|---|---|
| `ws` | 3 | Window size (px). Larger = fewer false tops, but merges close crowns |
| `method` | `"median"` | `median`, `mean`, `gaussian`, `maximum`. Median keeps crown edges sharp |

**Tree tops** — `detect_tops`, `merge_tops`, `screen_tops`

| Parameter | Default | Description |
|---|---|---|
| `hmin` | 2 | Minimum height (m) for a pixel to be a candidate |
| `ws` | 3 | Local-maximum window (px) = minimum spacing between tops |
| `distance` | 5 | Merge radius (px). Grouping is transitive — keep below crown diameter |

**Growing** — `grow_crowns`

| Parameter | Default | Description |
|---|---|---|
| `variance_thresh` | 2 | Max height variance within a crown. **The main lever** |
| `mask_thresh` | 0 | Minimum CHM height treated as canopy (m) |
| `morpho_radius` | 0 | Radius for opening/closing the mask. 0 = off |
| `alpha`, `beta`, `gamma` | 1, 0.5, 0.1 | Edge weights: mean diff, σ diff, inverse border length |
| `anneal_lambda` | 1 | Per-iteration tightening of the threshold. 1 = constant |
| `max_iters` | 200 | Cap on grow iterations per seed |
| `conflict_rule` | `"height"` | Who wins contested canopy — see below |
| `protect_seeds` | `FALSE` | If TRUE, no tree is absorbed; every top yields a crown |
| `retry_rejected` | `FALSE` | Reconsider regions rejected earlier in the same grow |

### Conflict rules

Two crowns can claim the same region. `conflict_rule` decides who gets it:

- **`"height"`** (default) — the taller tree wins. Dominant trees overtop their neighbours, so ambiguous canopy goes to the taller crown.
- **`"distance"`** — the nearest seed wins, by distance from the region centroid. Classic ITC behaviour; splits rather than merges.
- **`"similarity"`** — the tree whose apex height best matches the region's mean height wins.

Ties resolve to the lower crown id, so output is fully reproducible. `attr(crowns, "n_contested")` reports how many regions were disputed.

`protect_seeds = TRUE` disables merging entirely — every input top keeps its own crown. Use it when the tops are trusted, e.g. field-measured.

## Agreement with pyHRG

Two implementations of one method are only useful if they agree, so this is measured, not asserted.

| Test | Result |
|---|---|
| Watershed vs `scikit-image`, 200 random scenes incl. plateaus | **0 of 422,687 px differ** |
| Full pipeline vs pyHRG, 60 random scenes, all 3 conflict rules | **0 of 120,977 px differ** |
| Crown counts, contested counts | identical throughout |

Reproduce it yourself:

```bash
pip install pyhrg
python3 tools/generate_pyhrg_reference.py
Rscript tools/cross_validate_against_pyhrg.R
```

Two details make that agreement possible, and both are easy to get wrong in a port:

- **R stores matrices by column, Python by row.** The neighbour visiting order decides how plateaus resolve, so the C code is handed `t(matrix)`, whose column-major buffer is byte-identical to a row-major NumPy array. Skip that and the watershed silently transposes its tie-breaking.
- **Floating-point addition is not associative.** The pixel pass traverses in the same order as pyHRG, so the same heights are summed in the same sequence and the variances match to the last bit.

## Notes on behaviour

**The crown count can be lower than the tree-top count.** That is the point: merged trees leave gaps in the id sequence. Use `protect_seeds = TRUE` if you need one crown per top.

**`retry_rejected` only bites in the middle.** A rejected region can become admissible later, because absorbing a large homogeneous neighbour can *lower* a crown's variance. On a 1105-tree synthetic scene this changed nothing at `variance_thresh` 2 or 8, and nothing at 120, but changed ~17% of the raster at 20. The direction is not predictable.

**Where the C is, and why.** Two routines are compiled: the watershed, which is a sequential priority flood that no amount of vectorising fixes, and the pixel pass, which touches every cell. Everything else — the graph, the growing, the arbitration — is ordinary R, because it runs once per *tree*, not once per pixel, and that is the part you might actually want to read or change.

## Relation to PyCrown

The pipeline shape — smooth the CHM, find tops as local maxima, delineate crowns — comes from [PyCrown](https://github.com/manaakiwhenua/pycrown) (Zörner et al. 2018), via pyHRG. **The Dalponte & Coomes delineation that PyCrown re-implements is not part of rHRG**; hierarchical region growing is the only method here. If you want Dalponte in R, [lidR](https://github.com/r-lidar/lidR) and [itcSegment](https://cran.r-project.org/package=itcSegment) do it well and there is no reason to duplicate them.

## Testing

```r
remotes::install_github("igorpawelec/rhrg")
# or, from a clone:
#   R CMD build rhrg && R CMD check rhrg_0.1.0.tar.gz
```

## Citation

If you use rHRG in your research, please cite the software and the work it builds on:

> Pawelec, I. (2026). *rhrg: individual tree crown delineation from canopy height models by Hierarchical Region Growing* [R package]. https://github.com/igorpawelec/rhrg

**Upstream project.** The pipeline structure derives from PyCrown:

> Zörner, J., Dymond, J., Shepherd, J., Jolly, B. (2018). *PyCrown — Fast raster-based individual tree segmentation for LiDAR data.* Landcare Research NZ Ltd. https://doi.org/10.7931/M0SR-DN55
>
> Zörner, J., Dymond, J.R., Shepherd, J.D., Wiser, S.K., Bunting, P., Jolly, B. (2018). LiDAR-based regional inventory of tall trees — Wellington, New Zealand. *Forests* 9(11), 702. https://doi.org/10.3390/f9110702

**Methods used**

> Chan, T.F., Golub, G.H., LeVeque, R.J. (1983). Algorithms for computing the sample variance: analysis and recommendations. *The American Statistician* 37(3), 242–247. — the pairwise formula that merges two regions' statistics in O(1). Note this is *not* Welford's update, which adds one sample at a time; a crown absorbs a whole region at once.
>
> Beucher, S., Meyer, F. (1993). The morphological approach to segmentation: the watershed transformation. In: *Mathematical Morphology in Image Processing*, 433–481.
>
> Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest. *Photogrammetric Engineering & Remote Sensing* 70(5), 589–604. — local-maxima tree detection on a smoothed CHM.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

rHRG derives from PyCrown, which is published under GPLv3; the licence carries over.
