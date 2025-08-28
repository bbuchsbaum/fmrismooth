#' Fast joint bilateral using a high-dimensional lattice (splat/blur/slice)
#'
#' High-dimensional bilateral filtering using a permutohedral lattice backend.
#' Supports optional spatial guide and additional per-voxel guides.
#'
#' @param vol 3D volume array or `NeuroVol`.
#' @param sigma_sp Spatial sigma (voxels).
#' @param sigma_r Range sigma(s) for guide intensity dimensions. Scalar or
#'   vector matching `guide` + `guides` count.
#' @param guide Optional 3D spatial guide.
#' @param guides Optional list of additional 3D guides (e.g., probabilities).
#' @param blur_iters Integer number of lattice blur iterations.
#' @param mask Optional 3D logical/0-1 mask.
#' @return Smoothed 3D volume, wrapped like `vol` when possible.
#' @export
fast_bilateral_lattice3d <- function(vol,
                                     sigma_sp = 2.5,
                                     sigma_r  = 1.0,
                                     guide = NULL,
                                     guides = NULL,
                                     blur_iters = 2L,
                                     mask = NULL) {
  X <- nv_as_array(vol); X <- assert_numeric_array(X, 3L)
  if (!is.null(guide))  guide  <- .align_guide3d(guide, vol, kind = 'intensity')
  if (!is.null(guides)) guides <- .align_guides_list(guides, vol, kind = 'prob')
  G1 <- if (is.null(guide)) NULL else assert_numeric_array(nv_as_array(guide), 3L)
  GL <- if (is.null(guides)) NULL else lapply(guides, function(g) assert_numeric_array(nv_as_array(g), 3L))
  m  <- nv_check_mask(if (is.null(mask)) NULL else .align_mask_3d(mask, vol), dim(X))
  sr <- as.numeric(sigma_r); if (length(sr)==0) sr <- 1.0
  if (!is.null(GL) && length(sr) < length(GL) + (!is.null(G1))) {
    sr <- rep_len(sr, length(GL) + (!is.null(G1)))
  }
  Y <- .Call("_fmrismooth_bilateral_lattice3d_cpp", X, m,
             as.numeric(sigma_sp), as.numeric(sr),
             if (is.null(G1)) NULL else G1,
             if (is.null(GL)) NULL else GL,
             as.integer(blur_iters),
             PACKAGE = "fmrismooth")
  nv_wrap_like(vol, Y)
}

#' 4D joint bilateral lattice with weak temporal coupling (design-aware)
#'
#' Applies joint bilateral filtering across space and time. Temporal coupling
#' is modeled via an extra feature; an optional design regressor can be added
#' to the lattice features.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param sigma_sp Spatial sigma (voxels).
#' @param sigma_t Temporal sigma (frames).
#' @param sigma_r Range sigma(s) for guide intensity dimensions. Scalar or
#'   vector matching `guide_spatial` + length(`guides`).
#' @param guide_spatial Optional 3D spatial guide.
#' @param guides Optional list of 3D guides (probability maps, etc.).
#' @param design Optional numeric vector of length T (frames) for design-aware features.
#' @param sigma_d Sigma for design regressor feature.
#' @param blur_iters Integer number of lattice blur iterations.
#' @param mask Optional 3D logical/0-1 mask.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
fast_bilateral_lattice4d <- function(vec,
                                     sigma_sp = 2.5,
                                     sigma_t  = 0.5,
                                     sigma_r  = 1.0,
                                     guide_spatial = NULL,
                                     guides = NULL,
                                     design = NULL,
                                     sigma_d = 1.0,
                                     blur_iters = 1L,
                                     mask = NULL) {
  # Basic parameter validation
  if (!is.numeric(sigma_sp) || length(sigma_sp)!=1L || !is.finite(sigma_sp)) stop("sigma_sp must be a finite numeric scalar")
  if (!is.numeric(sigma_t)  || length(sigma_t)!=1L  || !is.finite(sigma_t))  stop("sigma_t must be a finite numeric scalar")
  if (!is.numeric(sigma_r)  || length(sigma_r)<1L)  stop("sigma_r must be numeric (>=1 value)")
  if (!is.numeric(sigma_d)  || length(sigma_d)!=1L || !is.finite(sigma_d))  stop("sigma_d must be a finite numeric scalar")
  if (!is.numeric(blur_iters) || length(blur_iters)!=1L) stop("blur_iters must be a scalar integer-like")
  X <- nv_as_array(vec); X <- assert_numeric_array(X, 4L)
  d <- dim(X)
  if (!is.null(guide_spatial)) guide_spatial <- .align_guide3d(guide_spatial, vec, kind='intensity')
  if (!is.null(guides))        guides        <- .align_guides_list(guides, vec, kind='prob')
  G1 <- if (is.null(guide_spatial)) NULL else assert_numeric_array(nv_as_array(guide_spatial), 3L)
  GL <- if (is.null(guides)) NULL else lapply(guides, function(g) assert_numeric_array(nv_as_array(g), 3L))
  m  <- nv_check_mask(if (is.null(mask)) NULL else .align_mask_3d(mask, vec), d)
  sr <- as.numeric(sigma_r); if (length(sr)==0) sr <- 1.0
  if (!is.null(GL) && length(sr) < length(GL) + (!is.null(G1))) {
    sr <- rep_len(sr, length(GL) + (!is.null(G1)))
  }
  if (!is.null(design)) {
    if (is.matrix(design) || is.data.frame(design)) {
      design <- as.matrix(design)
      if (nrow(design) != d[4]) stop("design matrix must have nrow = T (frames).")
    } else if (is.list(design)) {
      design <- lapply(design, as.numeric)
      if (any(vapply(design, length, 1L) != d[4])) stop("each design vector must have length T (frames).")
    } else {
      design <- as.numeric(design)
      if (length(design) != d[4]) stop("design must have length T (frames).")
    }
  }
  Y <- .Call("_fmrismooth_bilateral_lattice4d_cpp", X, m,
             as.numeric(sigma_sp), as.numeric(sigma_t), as.numeric(sr),
             if (is.null(G1)) NULL else G1,
             if (is.null(GL)) NULL else GL,
             if (is.null(design)) NULL else design,
             as.numeric(sigma_d),
             as.integer(blur_iters),
             PACKAGE = "fmrismooth")
  nv_wrap_like(vec, Y)
}
