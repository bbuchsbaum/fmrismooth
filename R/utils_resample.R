resample_to_target_space <- function(source, target, interpolation = 1L) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) return(NULL)
  sp_tgt <- tryCatch(neuroim2::space(target), error = function(e) NULL)
  if (is.null(sp_tgt)) return(NULL)
  if (inherits(source, "NeuroVol")) {
    rs <- tryCatch(neuroim2::resample(source, sp_tgt, interpolation = as.integer(interpolation)),
                   error = function(e) NULL)
    if (is.null(rs)) return(NULL)
    return(as.array(rs))
  }
  NULL
}

align_guide_3d <- function(guide, target, type = c("continuous","label"), interp = NULL) {
  type <- match.arg(type)
  if (is.null(interp)) interp <- if (type == "continuous") 1L else 0L
  ga <- if (inherits(guide, "NeuroVol")) as.array(guide) else nv_as_array(guide)
  td <- dim(nv_as_array(target))[1:3]
  if (length(dim(ga)) != 3L) stop("guide must be 3D")
  if (all(dim(ga) == td)) return(ga)
  if (inherits(guide, "NeuroVol")) {
    rs <- resample_to_target_space(guide, target, interpolation = interp)
    if (!is.null(rs) && all(dim(rs) == td)) return(rs)
  }
  stop("guide grid != target grid; supply NeuroVol with space or pre-align")
}

align_mask_3d <- function(mask, target, interp = 0L) {
  td <- dim(nv_as_array(target))[1:3]
  if (is.null(mask)) return(array(TRUE, dim = td))
  ma <- if (inherits(mask, "NeuroVol")) as.array(mask) else nv_as_array(mask)
  if (length(dim(ma)) > 3L) ma <- apply(ma, 1:3, function(v) any(v != 0))
  if (all(dim(ma) == td)) return(ma != 0)
  if (inherits(mask, "NeuroVol")) {
    rs <- resample_to_target_space(mask, target, interpolation = as.integer(interp))
    if (!is.null(rs) && all(dim(rs) == td)) return(rs != 0)
  }
  stop("mask grid != target grid; supply NeuroVol with space or pre-align")
}

# simple robust magnitude-MRI noise estimate from temporal differences (global version)
.estimate_sigma_rician_global <- function(vec, mask) {
  X <- nv_as_array(vec); stopifnot(length(dim(X)) == 4L)
  m <- if (is.null(mask)) array(TRUE, dim(X)[1:3]) else mask
  # std of temporal differences / sqrt(2)
  diffs <- X[,,, 2:dim(X)[4], drop=FALSE] - X[,,, 1:(dim(X)[4]-1), drop=FALSE]
  M4 <- array(rep(m, dim(diffs)[4]), dim = dim(diffs))
  s <- stats::mad(as.numeric(diffs[M4]), constant = 1.4826, na.rm = TRUE) / sqrt(2)
  s
}
