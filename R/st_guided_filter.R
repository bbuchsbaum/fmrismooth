#' 3D edge-preserving guided filter
#' @param vol 3D array or neuroim2::NeuroVol
#' @param radius integer window radius (voxels)
#' @param eps regularization (variance floor) in signal units^2
#' @param guide 3D array or NeuroVol; if NULL uses \code{vol}
#' @param mask 3D logical/0-1 mask; NULL => all voxels
#' @return same type as input (NeuroVol if given), otherwise 3D array
#' @export
st_guided_filter3d <- function(vol, radius = 4L, eps = 0.7^2,
                               guide = NULL, mask = NULL) {
  x  <- nv_as_array(vol); x <- assert_numeric_array(x, 3L)
  if (!is.null(guide)) guide <- .align_guide3d(guide, vol, kind='intensity')
  g  <- if (is.null(guide)) x else assert_numeric_array(nv_as_array(guide), 3L)
  m  <- nv_check_mask(if (is.null(mask)) NULL else .align_mask_3d(mask, vol), dim(x))
  stopifnot(is.finite(eps), eps >= 0, radius >= 0)
  out <- .Call("_fmrismooth_guided3d_core",
               x, g, m, as.integer(radius), as.numeric(eps),
               PACKAGE = "fmrismooth")
  nv_wrap_like(vol, out)
}

#' 4D separable spatiotemporal guided filter
#' Spatial guided filtering is applied per-frame (3D). Optionally, a temporal
#' guided step (1D per voxel) is applied with motion-aware frame weights.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param radius_sp Spatial filter radius (voxels).
#' @param eps_sp Spatial regularization (variance floor).
#' @param radius_t Temporal filter radius (frames).
#' @param eps_t Temporal regularization.
#' @param guide_spatial Optional 3D spatial guide.
#' @param guide_temporal Optional length-T temporal guide.
#' @param motion_params Optional T×(3 or 6) motion regressors; used to weight temporal guidance.
#' @param sigma_m Motion weight scale (frames).
#' @param mask Optional 3D logical/0-1 mask.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
st_guided_filter4d <- function(vec,
                               radius_sp = 4L, eps_sp = 0.7^2,
                               radius_t = 1L, eps_t = 0.3^2,
                               guide_spatial = NULL,
                               guide_temporal = NULL,
                               motion_params = NULL, sigma_m = 0.5,
                               mask = NULL) {
  X <- nv_as_array(vec); X <- assert_numeric_array(X, 4L)
  d <- dim(X); nt <- d[4]
  if (!is.null(guide_spatial)) guide_spatial <- .align_guide3d(guide_spatial, vec, kind='intensity')
  m <- nv_check_mask(if (is.null(mask)) NULL else .align_mask_3d(mask, vec), d)

  # 1) Spatial guided, per-frame
  Gs <- if (is.null(guide_spatial)) NULL else assert_numeric_array(nv_as_array(guide_spatial), 3L)
  for (t in seq_len(nt)) {
    X[,,,t] <- .Call("_fmrismooth_guided3d_core",
                     X[,,,t], if (is.null(Gs)) X[,,,t] else Gs, m,
                     as.integer(radius_sp), as.numeric(eps_sp),
                     PACKAGE = "fmrismooth")
  }

  # 2) Temporal guided, per-voxel (optional)
  if (radius_t > 0L) {
    if (!is.null(motion_params)) {
      mp <- as.matrix(motion_params)
      if (!ncol(mp) %in% c(3,6)) stop("motion_params should have 6 (or 3) cols.")
      if (nrow(mp) != nt) stop("motion_params rows must equal number of frames.")
      dmp <- sqrt(rowSums(rbind(rep(0, ncol(mp)), diff(mp))^2))
      w   <- exp(-(dmp^2) / (2 * sigma_m^2))
    } else {
      w <- rep(1.0, nt)
    }
    if (is.null(guide_temporal)) {
      X <- .Call("_fmrismooth_temporal_guided1d_core",
                 X, NULL, as.numeric(w), m,
                 as.integer(radius_t), as.numeric(eps_t),
                 PACKAGE = "fmrismooth")
    } else {
      gt <- as.numeric(guide_temporal)
      if (length(gt) != nt) stop("guide_temporal must have length T (number of frames).")
      X <- .Call("_fmrismooth_temporal_guided1d_core",
                 X, gt, as.numeric(w), m,
                 as.integer(radius_t), as.numeric(eps_t),
                 PACKAGE = "fmrismooth")
    }
  }
  nv_wrap_like(vec, X)
}
