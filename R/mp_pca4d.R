#' MP-PCA denoising (local patches × time) for 4D fMRI
#'
#' Performs local PCA denoising in overlapping 3D spatial patches across a
#' temporal window.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param patch Integer length-3 patch size (voxels).
#' @param tw Temporal window size (frames).
#' @param stride Integer length-4 stride for x,y,z,t.
#' @param sigma_mode Noise mode: 'patch' (estimate per-patch), 'global' (one value), or 'fixed'.
#' @param sigma_value Noise sigma when `sigma_mode='fixed'`.
#' @param mask Optional 3D logical/0-1 mask.
#' @return Denoised 4D data, wrapped like `vec` when possible.
#' @export
mp_pca4d <- function(vec,
                     patch = c(5L,5L,5L),
                     tw = 32L,
                     stride = c(3L,3L,3L,8L),
                     sigma_mode = c('patch','global','fixed'),
                     sigma_value = NULL,
                     mask = NULL) {
  X <- nv_as_array(vec); stopifnot(length(dim(X))==4L)
  d <- dim(X); m <- nv_check_mask(mask, d)
  patch <- as.integer(patch); if (length(patch)!=3L) stop("patch must be length-3")
  tw <- as.integer(tw); stride <- as.integer(stride); if (length(stride)!=4L) stop("stride must be length-4")
  sigma_mode <- match.arg(sigma_mode)
  if (sigma_mode == 'patch') {
    sig <- NA_real_
  } else if (sigma_mode == 'global') {
    sig <- estimate_sigma_rician(X, mask = m)
  } else {
    if (is.null(sigma_value)) stop("sigma_value must be provided when sigma_mode='fixed'")
    sig <- as.numeric(sigma_value)
  }
  Y <- .Call("_fmrismooth_mppca4d_core",
             X, m, patch, as.integer(tw), stride, sig,
             PACKAGE = "fmrismooth")
  nv_wrap_like(vec, Y)
}
