#' Robust noise sigma estimate for Rician magnitude fMRI (4D only)
#'
#' Estimates the global noise standard deviation from temporal differences,
#' assuming Rician magnitude noise that is approximately Gaussian at moderate SNR.
#'
#' @param vec 4D fMRI data as a numeric array or `neuroim2::NeuroVec`.
#' @param mask Optional 3D logical/0-1 mask (array or `NeuroVol`). If `NULL`,
#'   all voxels are used.
#' @return A single numeric value: estimated noise sigma.
#' @examples
#' set.seed(1)
#' x <- array(100 + rnorm(10*10*10*5, sd=2), c(10,10,10,5))
#' estimate_sigma_rician(x)
#' @export
estimate_sigma_rician <- function(vec, mask = NULL) {
  X <- nv_as_array(vec); stopifnot(length(dim(X))==4L)
  m <- nv_check_mask(mask, dim(X))
  nt <- dim(X)[4]
  if (nt < 3L) stop("Need at least 3 frames to estimate sigma.")
  dx <- X[,,,2:nt] - X[,,,1:(nt-1)]
  mad_scale <- apply(dx, 1:3, function(ts) stats::mad(ts, center = stats::median(ts), constant = 1.0))
  sd_est <- (1.4826 * mad_scale) / sqrt(2.0)
  sd_est_vec <- as.numeric(sd_est[m])
  stats::median(sd_est_vec, na.rm = TRUE)
}

#' VST forward for Rician magnitude
#'
#' Applies a variance-stabilizing transform to Rician magnitude data.
#'
#' @param x Numeric array (3D/4D) of magnitudes.
#' @param sigma Noise standard deviation.
#' @return Array of the same shape as `x`.
#' @export
vst_forward_rician <- function(x, sigma) {
  xc <- pmax(x, 0)
  sqrt(pmax(xc * xc - 2 * sigma^2, 0))
}

#' VST inverse for Rician magnitude
#'
#' Inverts the Rician VST.
#'
#' @param z Transformed array (output of `vst_forward_rician`).
#' @param sigma Noise standard deviation.
#' @return Array of the same shape as `z`.
#' @export
vst_inverse_rician <- function(z, sigma) {
  zc <- pmax(z, 0)
  sqrt(pmax(zc * zc + 2 * sigma^2, 0))
}

#' Wrap a denoiser with a Rician VST (variance-stabilizing transform)
#'
#' Converts magnitude data to approximately homoscedastic Gaussian via VST,
#' applies `denoise_fun`, then inverts the VST.
#'
#' @param x 3D/4D magnitude array or neuroim2 volume/vec.
#' @param sigma Optional noise sigma. If `NULL` and a 4D input is provided,
#'   sigma is estimated with `estimate_sigma_rician`.
#' @param denoise_fun Function that takes an array and returns an array of the
#'   same dimensions.
#' @param ... Additional arguments forwarded to `denoise_fun`.
#' @return The denoised data, same shape/type as `x` (wrapped via `nv_wrap_like`).
#' @export
vst_wrap <- function(x, sigma = NULL, denoise_fun, ...) {
  arr <- nv_as_array(x)
  if (is.null(sigma)) {
    if (length(dim(arr)) != 4L)
      stop("sigma must be supplied for 3D; can be estimated for 4D.")
    sigma <- estimate_sigma_rician(arr)
  }
  z <- vst_forward_rician(arr, sigma)
  z_dn <- denoise_fun(z, ...)
  # Ensure the denoiser preserves shape
  if (!is.array(z_dn) || !identical(dim(z_dn), dim(arr))) {
    stop("denoise_fun must return an array with the same dimensions as the input.")
  }
  out <- vst_inverse_rician(z_dn, sigma)
  nv_wrap_like(x, out)
}
