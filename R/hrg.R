#' @useDynLib rhrg, .registration = TRUE
NULL

#' Conflict-resolution rules
#'
#' The rules [grow_crowns()] accepts for canopy claimed by more than one
#' crown: `"height"`, `"distance"`, `"similarity"`.
#' @export
CONFLICT_RULES <- c("height", "distance", "similarity")

# Combine two sets of (count, mean, population variance) in O(1) using the
# pairwise formula of Chan, Golub & LeVeque (1983), which merges subsamples
# of arbitrary size. Welford's update is a different algorithm: it adds one
# sample at a time, whereas a crown absorbs a whole region at once.
.merge_stats <- function(nA, meanA, varA, nB, meanB, varB) {
  N <- nA + nB
  if (N == 0) return(c(0, 0, 0))
  mean_new <- (nA * meanA + nB * meanB) / N
  var_new <- (nA * (varA + (meanA - mean_new)^2) +
              nB * (varB + (meanB - mean_new)^2)) / N
  c(N, mean_new, var_new)
}

#' Grow tree crowns from a CHM and a set of tree tops
#'
#' The algorithm in three stages. A marker-based watershed on the inverted
#' CHM, seeded with the tops, yields exactly one region per top, so the
#' regions *are* the detected trees. Neighbouring regions are linked in a
#' weighted graph. Each tree then grows greedily, absorbing neighbours while
#' the combined height variance stays below `variance_thresh`; when two trees
#' absorb each other they are one crown, which is how over-detected tops get
#' corrected.
#'
#' Growing runs independently per seed, so two crowns can claim the same
#' region. `conflict_rule` settles those claims on the data rather than on
#' iteration order, which keeps the result reproducible whatever order the
#' tops arrive in.
#'
#' @param chm Numeric matrix, ideally already smoothed.
#' @param tops Two-column matrix of integer (row, col) tree top positions,
#'   1-based, as returned by [as_pixels()].
#' @param variance_thresh Numeric. Maximum height variance within a crown.
#'   The main control: raise it and neighbouring trees merge more readily.
#' @param mask_thresh Numeric. Minimum CHM height treated as canopy.
#' @param morpho_radius Integer. Radius for opening then closing the mask,
#'   to drop speckle. 0 disables it.
#' @param alpha,beta,gamma Numeric. Edge weights for the mean difference,
#'   the standard-deviation difference and the inverse shared border length.
#' @param anneal_lambda Numeric. Per-iteration tightening of the variance
#'   threshold. 1 keeps it constant.
#' @param max_iters Integer. Cap on grow iterations per seed.
#' @param conflict_rule One of [CONFLICT_RULES]. `"height"` gives contested
#'   canopy to the taller tree, `"distance"` to the nearest, `"similarity"`
#'   to the tree whose apex height best matches the region. Ties resolve to
#'   the lower crown id, so output is reproducible.
#' @param protect_seeds Logical. If `TRUE` no tree is ever absorbed and every
#'   top yields a crown, disabling merging.
#' @param retry_rejected Logical. Reconsider regions rejected earlier in the
#'   same grow. Only bites at intermediate `variance_thresh`, where crowns
#'   actively merge.
#'
#' @return An integer matrix of crown ids, 0 for background, same dimensions
#'   as `chm`. Ids follow the order of `tops`; absorbed trees leave gaps, so
#'   the crown count can be lower than the number of tops. The number of
#'   contested regions is attached as attribute `"n_contested"`.
#'
#' @references
#' Chan, T.F., Golub, G.H., LeVeque, R.J. (1983). Algorithms for computing
#' the sample variance. *The American Statistician* 37(3), 242-247.
#'
#' Beucher, S., Meyer, F. (1993). The morphological approach to
#' segmentation: the watershed transformation.
#'
#' @examples
#' yy <- row(matrix(0, 40, 40)); xx <- col(matrix(0, 40, 40))
#' chm <- 20 * exp(-((yy - 20)^2 + (xx - 20)^2) / 60)
#' grow_crowns(chm, matrix(c(20L, 20L), ncol = 2), mask_thresh = 1)
#' @export
grow_crowns <- function(chm, tops,
                        variance_thresh = 2,
                        mask_thresh = 0,
                        morpho_radius = 0L,
                        alpha = 1, beta = 0.5, gamma = 0.1,
                        anneal_lambda = 1,
                        max_iters = 200L,
                        conflict_rule = "height",
                        protect_seeds = FALSE,
                        retry_rejected = FALSE) {

  chm <- as.matrix(chm)
  if (length(dim(chm)) != 2L) stop("chm must be a 2-D matrix", call. = FALSE)
  if (!conflict_rule %in% CONFLICT_RULES)
    stop("conflict_rule must be one of ", paste(CONFLICT_RULES, collapse = ", "),
         call. = FALSE)

  nr <- nrow(chm); nc <- ncol(chm)
  tops <- matrix(as.integer(tops), ncol = 2)
  if (nrow(tops) > 0L &&
      any(tops[, 1] < 1L | tops[, 1] > nr | tops[, 2] < 1L | tops[, 2] > nc))
    stop("a tree top lies outside the CHM (", nr, "x", nc, ")", call. = FALSE)

  n_seeds <- nrow(tops)
  markers <- matrix(0L, nr, nc)
  if (n_seeds > 0L) markers[tops] <- seq_len(n_seeds)
  mask <- .prepare_mask(chm, mask_thresh, morpho_radius)

  labels <- .watershed(-chm, markers, mask)
  if (n_seeds == 0L) {
    out <- matrix(0L, nr, nc); attr(out, "n_contested") <- 0L; return(out)
  }

  g <- .rag(labels, chm, n_seeds)
  csr <- .csr(g, n_seeds, alpha, beta, gamma)

  members <- vector("list", n_seeds)
  for (s in seq_len(n_seeds))
    members[[s]] <- .grow_one(s, csr, g, variance_thresh, max_iters,
                              anneal_lambda, retry_rejected)

  seed_h <- chm[tops]
  res <- .resolve_conflicts(members, tops, seed_h, g, conflict_rule, protect_seeds)

  lut <- integer(max(labels) + 1L)
  lut[as.integer(names(res$assignment)) + 1L] <- res$assignment
  out <- matrix(lut[labels + 1L], nr, nc)
  attr(out, "n_contested") <- res$n_contested
  out
}

.watershed <- function(img, markers, mask) {
  nr <- nrow(img); nc <- ncol(img)
  # t() so the column-major buffer matches the row-major layout the C code
  # and pyHRG both assume; without it the neighbour order would transpose.
  out <- .Call(C_watershed, as.double(t(img)), as.integer(t(markers)),
               as.integer(t(mask)), nr, nc)
  matrix(out, nr, nc, byrow = TRUE)
}

.rag <- function(labels, chm, nreg) {
  .Call(C_rag, as.integer(t(labels)), as.double(t(chm)),
        nrow(labels), ncol(labels), as.integer(nreg))
}

.prepare_mask <- function(chm, thresh, radius) {
  m <- chm > thresh
  if (radius > 0L) {
    se <- .disk(radius)
    m <- .morph(m, se, min)   # erosion then dilation = opening
    m <- .morph(m, se, max)
    m <- .morph(m, se, max)   # then closing, to refill crown interiors
    m <- .morph(m, se, min)
  }
  matrix(as.integer(m), nrow(chm), ncol(chm))
}

.disk <- function(r) {
  g <- -r:r
  outer(g, g, function(a, b) a^2 + b^2 <= r^2)
}

.morph <- function(m, se, fun) {
  r <- (nrow(se) - 1L) %/% 2L
  nr <- nrow(m); nc <- ncol(m)
  pm <- matrix(if (identical(fun, min)) TRUE else FALSE, nr + 2L*r, nc + 2L*r)
  pm[seq_len(nr) + r, seq_len(nc) + r] <- m
  out <- matrix(if (identical(fun, min)) TRUE else FALSE, nr, nc)
  first <- TRUE
  for (i in seq_len(nrow(se))) for (j in seq_len(ncol(se))) {
    if (!se[i, j]) next
    w <- pm[seq_len(nr) + i - 1L, seq_len(nc) + j - 1L, drop = FALSE]
    out <- if (first) w else if (identical(fun, min)) out & w else out | w
    first <- FALSE
  }
  out
}

# CSR adjacency with weighted edges:
#   w(a,b) = alpha*|dmu| + beta*|dsigma| + gamma/(shared border + 1)
.csr <- function(g, nreg, alpha, beta, gamma) {
  ea <- g$ea; eb <- g$eb
  if (length(ea) == 0L)
    return(list(row_ptr = integer(nreg + 2L), col_idx = integer(0),
                weights = numeric(0)))
  key <- ea * (nreg + 1) + eb
  ub <- !duplicated(key)
  border <- as.vector(table(key))[match(key[ub], sort(unique(key)))]
  ua <- ea[ub]; ubb <- eb[ub]

  # undirected: each pair appears once per direction
  from <- c(ua, ubb); to <- c(ubb, ua); bl <- c(border, border)
  o <- order(from, to)
  from <- from[o]; to <- to[o]; bl <- bl[o]

  sdv <- sqrt(pmax(g$var, 0))
  w <- alpha * abs(g$mean[from + 1L] - g$mean[to + 1L]) +
       beta  * abs(sdv[from + 1L] - sdv[to + 1L]) +
       gamma * (1 / (bl + 1))

  # CSR offsets: row_ptr[k+1] is where node k's neighbours start. Nodes are
  # 1-based, so with R's 1-based vectors row_ptr[1] and row_ptr[2] are both 0.
  row_ptr <- integer(nreg + 2L)
  cnt <- tabulate(from, nbins = nreg)
  row_ptr[3:(nreg + 2L)] <- cumsum(cnt)
  list(row_ptr = row_ptr, col_idx = to, weights = w)
}

# Greedy grow for one seed. Candidates come off a min-heap keyed by
# (weight, id), so ties resolve by the lower id and the walk is reproducible.
.grow_one <- function(seed, csr, g, v_thresh, max_iters, anneal_lambda,
                      retry_rejected) {
  cur <- c(g$n[seed + 1L], g$mean[seed + 1L], g$var[seed + 1L])
  members <- seed
  in_members <- rep(FALSE, length(g$n))
  in_members[seed + 1L] <- TRUE
  rejected <- rep(FALSE, length(g$n))
  vt <- v_thresh

  hw <- numeric(0); hi <- integer(0)
  push <- function(w, id) { hw <<- c(hw, w); hi <<- c(hi, id) }
  pop <- function() {
    k <- which.min(hw + 1e-12 * hi)     # (weight, id) ordering
    o <- order(hw, hi)[1]
    id <- hi[o]; hw <<- hw[-o]; hi <<- hi[-o]
    id
  }
  rng <- (csr$row_ptr[seed + 1L] + 1L):csr$row_ptr[seed + 2L]
  if (csr$row_ptr[seed + 2L] > csr$row_ptr[seed + 1L])
    for (k in rng) push(csr$weights[k], csr$col_idx[k])

  for (it in seq_len(max_iters)) {
    if (!length(hi)) break
    if (anneal_lambda < 1 && it > 1L) vt <- vt * anneal_lambda
    cand <- NA_integer_
    while (length(hi)) {
      id <- pop()
      if (in_members[id + 1L]) next
      if (rejected[id + 1L] && !retry_rejected) next
      cand <- id; break
    }
    if (is.na(cand)) break

    tst <- .merge_stats(cur[1], cur[2], cur[3],
                        g$n[cand + 1L], g$mean[cand + 1L], g$var[cand + 1L])
    if (tst[3] <= vt) {
      members <- c(members, cand)
      in_members[cand + 1L] <- TRUE
      rejected[cand + 1L] <- FALSE
      cur <- tst
      if (csr$row_ptr[cand + 2L] > csr$row_ptr[cand + 1L])
        for (k in (csr$row_ptr[cand + 1L] + 1L):csr$row_ptr[cand + 2L]) {
          nb <- csr$col_idx[k]
          if (!in_members[nb + 1L]) push(csr$weights[k], nb)
        }
    } else {
      rejected[cand + 1L] <- TRUE
    }
  }
  members
}

# Assign every region to exactly one crown, deciding contested ones on the
# data. See the package README for why write order is not good enough.
.resolve_conflicts <- function(members, tops, seed_h, g, rule, protect_seeds) {
  claims <- list()
  for (cid in seq_along(members))
    for (m in members[[cid]]) {
      k <- as.character(m)
      claims[[k]] <- c(claims[[k]], cid)
    }
  assignment <- integer(length(claims))
  names(assignment) <- names(claims)
  n_contested <- 0L

  for (k in names(claims)) {
    cl <- claims[[k]]
    if (length(cl) == 1L) { assignment[k] <- cl; next }
    n_contested <- n_contested + 1L
    ws <- as.integer(k)
    if (protect_seeds && ws %in% cl) { assignment[k] <- ws; next }
    score <- switch(rule,
      height     = -seed_h[cl],
      distance   = (g$cy[ws + 1L] - (tops[cl, 1] - 1))^2 +
                   (g$cx[ws + 1L] - (tops[cl, 2] - 1))^2,
      similarity = abs(seed_h[cl] - g$mean[ws + 1L]))
    assignment[k] <- cl[order(score, cl)[1]]     # ties -> lower crown id
  }
  list(assignment = assignment, n_contested = n_contested)
}
