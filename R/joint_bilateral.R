#' 3D joint bilateral (permutohedral backend; 'grid' alias supported)
#'
#' Joint bilateral filter using a permutohedral lattice backend (the `grid`
#' alias is supported and routed to the same implementation).
#'
#' @param vol 3D volume array or `NeuroVol`.
#' @param guide_spatial 3D guide image (array or `NeuroVol`).
#' @param sigma_s Spatial sigma (voxels).
#' @param sigma_r Range sigma for guide intensity.
#' @param range_bins Unused placeholder for grid backend compatibility.
#' @param passes Integer number of lattice blur passes.
#' @param mask Optional 3D mask.
#' @param interp_guide Guide interpolation mode for resampling.
#' @param interp_mask Mask interpolation mode for resampling.
#' @param backend One of 'grid' or 'permutohedral'.
#' @param guides Optional list of additional 3D guides.
#' @return Smoothed 3D data, wrapped like `vol` when possible.
#' @export
fast_bilateral_joint3d <- function(
  vol, guide_spatial,
  sigma_s, sigma_r, range_bins = NULL, passes = 1L,
  mask = NULL, interp_guide = 1L, interp_mask = 0L,
  backend = c("grid","permutohedral"),
  guides = NULL
) {
  backend <- match.arg(backend)
  X <- nv_as_array(vol); stopifnot(length(dim(X)) == 3L)
  G <- align_guide_3d(guide_spatial, vol, type = "continuous", interp = interp_guide)
  m <- align_mask_3d(mask, vol, interp = interp_mask)
  if (!all(dim(G) == dim(X))) stop("guide and vol must share 3D dims after alignment")
  # we route both backends to permutohedral (grid alias) in this build
  Y <- .Call('_fmrismooth_permutohedral_lattice3d_cpp', X, m,
             as.numeric(sigma_s), as.numeric(sigma_r),
             G,
             if (is.null(guides)) NULL else guides,
             as.integer(passes),
             PACKAGE = 'fmrismooth')
  nv_wrap_like(vol, Y)
}

#' 4D joint bilateral (permutohedral backend; 'grid' alias supported)
#'
#' Joint bilateral filter across space and time using a lattice backend.
#'
#' @param vec 4D array or `NeuroVec`.
#' @param guide_spatial 3D spatial guide.
#' @param sigma_s Spatial sigma (voxels).
#' @param sigma_t Temporal sigma (frames).
#' @param sigma_r Range sigma for guide intensity.
#' @param range_bins Unused placeholder for grid backend compatibility.
#' @param passes Integer number of lattice blur passes.
#' @param mask Optional 3D mask.
#' @param interp_guide Guide interpolation mode for resampling.
#' @param interp_mask Mask interpolation mode for resampling.
#' @param backend One of 'grid' or 'permutohedral'.
#' @param design Optional length-T design regressor to add as lattice feature.
#' @param sigma_d Sigma for the design feature.
#' @param guides Optional list of additional 3D guides.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
fast_bilateral_joint4d <- function(
  vec, guide_spatial, sigma_s, sigma_t, sigma_r,
  range_bins = NULL, passes = 1L, mask = NULL,
  interp_guide = 1L, interp_mask = 0L,
  backend = c("grid","permutohedral"),
  design = NULL, sigma_d = 1.0, guides = NULL
) {
  backend <- match.arg(backend)
  X <- nv_as_array(vec); stopifnot(length(dim(X)) == 4L)
  G <- align_guide_3d(guide_spatial, vec, type = "continuous", interp = interp_guide)
  m <- align_mask_3d(mask, vec, interp = interp_mask)
  if (!all(dim(G) == dim(X)[1:3])) stop("guide_spatial must match spatial dims after alignment")
  Y <- .Call('_fmrismooth_permutohedral_lattice4d_cpp', X, m,
             as.numeric(sigma_s), as.numeric(sigma_t), as.numeric(sigma_r),
             G,
             if (is.null(guides)) NULL else guides,
             if (is.null(design)) NULL else design,
             as.numeric(sigma_d),
             as.integer(passes),
             PACKAGE = 'fmrismooth')
  nv_wrap_like(vec, Y)
}
