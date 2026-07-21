# cran-comments

## Test environments

* Windows 11, R 4.5.2, Rtools45 (local), `R CMD check --as-cran`, with TinyTeX
  installed so that `checking PDF version of manual` runs rather than
  being skipped
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)

## R CMD check results

0 errors, 0 warnings, 1 note.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Igor Pawelec <igor.pawelec@urk.edu.pl>'
New submission
```

Expected for a first submission.

The remaining local notes are environmental rather than properties of the
package: `checking top-level files` reports that README.md and NEWS.md
cannot be checked without pandoc, and `checking for future file
timestamps` reports that it could not reach the network time service.
Neither appears on the CI runners, which install pandoc.

`checking PDF version of manual` passes. Note that the CI workflow runs
with `--no-manual` to avoid requiring a LaTeX install on five runners, so
that particular check is covered locally rather than in CI.

## Notes for the reviewer

The package compiles two C routines and uses no Rcpp: a marker-based
watershed, which is a sequential priority flood, and a single pixel pass
that accumulates per-region statistics. Everything else is plain R.
`stats` is the only import. `terra` is a Suggests, used only by
`read_chm()`, `write_crowns()` and `crowns_to_polygons()`; the algorithm
itself runs on ordinary numeric matrices.

The DESCRIPTION states that results match the pyHRG Python package pixel
for pixel on the shared synthetic test suite, and closely but not exactly
on canopy height models containing plateaus. That distinction is
deliberate and measured, not hedging: the two implementations break
plateau ties differently, the difference is 0.25% of watershed pixels on
the real scenes tested, and both figures are documented in the README
along with how to reproduce them.

`tools/` is excluded from the build via `.Rbuildignore`. It holds the
cross-check against pyHRG, which needs Python and so cannot run on a CRAN
machine, and a script that measures the method's sensitivity to those
plateau ties.
