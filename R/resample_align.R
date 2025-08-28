# neuroim2-backed resampling of guides/masks when grids differ

#' @keywords internal
.has_neuroim2 <- function() requireNamespace("neuroim2", quietly = TRUE)

#' @keywords internal
.is_neurovol <- function(x) inherits(x, "NeuroVol")

#' @keywords internal
.is_neurovec <- function(x) inherits(x, "NeuroVec")

#' @keywords internal
.get_space <- function(x) {
  if (.has_neuroim2() && (inherits(x, "NeuroVol") || inherits(x, "NeuroVec"))) {
    return(neuroim2::space(x))
  }
  NULL
}

#' @keywords internal
.align_3d_to_target_space <- function(obj, target, kind = c("intensity","mask","prob")) {
  kind <- match.arg(kind)
  # Coerce numeric vector to 3D array if it matches target's spatial size
  if (!is.array(obj) && is.numeric(obj)) {
    td <- dim(target); if (length(td) == 4L) td <- td[1:3]
    if (length(td) == 3L && length(obj) == prod(td)) {
      obj <- array(obj, dim = td)
    }
  }
  # If both are arrays, accept 3D obj against 3D or 4D target when spatial dims match
  if (is.array(obj) && is.array(target)) {
    od <- dim(obj)
    td <- dim(target)
    if (length(od) == 3L) {
      if (length(td) == 4L) td <- td[1:3]
      if (length(td) == 3L && all(od == td)) return(obj)
      stop("Guides/masks have different array dimensions and no spatial info; install 'neuroim2' and pass NeuroVol/NeuroVec to enable resampling.")
    }
  }
  # If target has spatial metadata and obj is NeuroVol, try resampling via neuroim2
  trg_space <- .get_space(target)
  if (.has_neuroim2() && .is_neurovol(obj) && !is.null(trg_space)) {
    interp <- switch(kind, mask=0L, prob=1L, intensity=1L, 1L)
    out <- tryCatch(neuroim2::resample(obj, trg_space, interpolation = interp),
                    error = function(e) stop("neuroim2::resample failed: ", conditionMessage(e)))
    return(as.array(out))
  }
  # Final permissive check if both are arrays and spatial dims match
  if (is.array(obj) && is.array(target)) {
    od <- dim(obj); td <- dim(target); if (length(td) == 4L) td <- td[1:3]
    if (length(od) == 3L && length(td) == 3L && all(od == td)) return(obj)
  }
  stop("Could not align guide/mask to target grid. Install 'neuroim2 (>= 0.8.1)' and pass NeuroVol/NeuroVec to enable resampling.")
}

#' @keywords internal
.align_mask_3d <- function(mask, target) {
  if (is.null(mask)) return(NULL)
  .align_3d_to_target_space(mask, target, kind = "mask")
}

#' @keywords internal
.align_guide3d <- function(guide, target, kind = c("intensity","prob")) {
  if (is.null(guide)) return(NULL)
  .align_3d_to_target_space(guide, target, kind = match.arg(kind))
}

#' @keywords internal
.align_guides_list <- function(guides, target, kind = c("intensity","prob")) {
  if (is.null(guides)) return(NULL)
  kind <- match.arg(kind)
  lapply(guides, function(g) .align_guide3d(g, target, kind = kind))
}
