#' fmrismooth_mppca_pipeline: MP-PCA denoising followed by bilateral filtering
#'
#' Denoises 4D fMRI with MP-PCA in overlapping patches, then applies a
#' joint-bilateral lattice filter with optional guides.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param t1 Optional 3D structural guide.
#' @param probs Optional list of 3D guides (e.g., tissue probabilities).
#' @param mask Optional 3D logical/0-1 mask.
#' @param design Optional per-frame design regressor.
#' @param motion Optional per-frame motion matrix used by some stages.
#' @param sigma_mode MP-PCA noise mode: 'patch', 'global', or 'fixed'.
#' @param sigma_value Noise sigma when `sigma_mode='fixed'`.
#' @param sigma_sp Spatial sigma (voxels) for joint bilateral.
#' @param sigma_t Temporal sigma (frames) for joint bilateral.
#' @param sigma_r Range sigma(s) for joint bilateral guidance.
#' @param sigma_d Sigma for design regressor feature.
#' @param lattice_blur_iters Number of lattice blur iterations.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
fmrismooth_mppca_pipeline <- function(vec,
                               t1 = NULL,
                               probs = NULL,
                               mask = NULL,
                               design = NULL,
                               motion = NULL,
                               sigma_mode = c('patch','global','fixed'),
                               sigma_value = NULL,
                               sigma_sp = 2.5,
                               sigma_t  = 0.5,
                               sigma_r  = 15,
                               sigma_d  = 1.0,
                               lattice_blur_iters = 1L) {
  X <- nv_as_array(vec); stopifnot(length(dim(X))==4L)
  d <- dim(X)
  if (!is.null(t1))    t1    <- .align_guide3d(t1, vec, kind='intensity')
  if (!is.null(probs)) probs <- .align_guides_list(probs, vec, kind='prob')
  m <- nv_check_mask(if (is.null(mask)) NULL else .align_mask_3d(mask, vec), d)

  Xden <- mp_pca4d(vec, sigma_mode = sigma_mode[1], sigma_value = sigma_value,
                   patch = c(5L,5L,5L), tw = min(32L, d[4]), stride = c(3L,3L,3L,8L), mask = m)

  out <- fast_bilateral_lattice4d(
    Xden,
    sigma_sp = sigma_sp, sigma_t = sigma_t, sigma_r = sigma_r,
    guide_spatial = t1, guides = probs, design = design, sigma_d = sigma_d,
    blur_iters = lattice_blur_iters, mask = m
  )
  nv_wrap_like(vec, out)
}
