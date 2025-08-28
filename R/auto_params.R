#' Recommend smoothing parameters from voxel size, TR, and noise
#'
#' Heuristically selects spatial/temporal smoothing and robust-TV weights based
#' on voxel spacing, repetition time, and a noise estimate.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param mask Optional 3D logical/0-1 mask.
#' @param tr Optional repetition time (seconds). Defaults to 2.
#' @param target_fwhm_mm Target spatial FWHM in millimeters.
#' @param sigma_mppca Optional pre-computed noise sigma. When `NULL`, noise is
#'   estimated with `estimate_sigma_rician`.
#' @return A list with fields `lambda_s`, `lambda_t`, `sigma_sp`, `sigma_t`, `sigma_r`.
#' @export
recommend_params <- function(vec, mask = NULL, tr = NULL, target_fwhm_mm = 5, sigma_mppca = NULL) {
  X <- nv_as_array(vec); stopifnot(length(dim(X))==4L)
  d <- dim(X)
  m <- tryCatch(align_mask_3d(mask, vec), error=function(e) array(TRUE, d[1:3]))
  vx <- c(3,3,3)
  if (requireNamespace("neuroim2", quietly = TRUE) && inherits(vec, "NeuroVec")) {
    sp <- tryCatch(neuroim2::space(vec), error = function(e) NULL)
    if (!is.null(sp)) vx <- tryCatch(neuroim2::spacing(sp), error = function(e) c(3,3,3))
  }
  tr <- if (is.null(tr)) 2.0 else as.numeric(tr)
  sig <- if (!is.null(sigma_mppca)) as.numeric(sigma_mppca) else estimate_sigma_rician(X, mask = m)
  # take median over all timepoints within the 3D mask
  M4 <- array(rep(m, d[4]), dim = d)
  mu  <- stats::median(as.numeric(X[M4]), na.rm = TRUE)
  sdr <- stats::mad(as.numeric(apply(X, 1:3, mean)[m]), constant = 1.4826, na.rm = TRUE)
  fwhm_vox <- target_fwhm_mm / vx[1]
  sigma_sp <- max(1.0, (fwhm_vox / 2.355))
  sigma_t  <- max(0.3, min(0.8, 0.5 * (tr / 2.0)))
  sigma_r  <- max(5, 1.5 * sdr)
  k_s <- 0.7; k_t <- 0.006
  lambda_s <- k_s * (mu / (sig + 1e-6)) * (vx[1]/mean(vx))
  lambda_t <- k_t * (mu / (sig + 1e-6)) * (2.0 / tr)
  lambda_s <- max(0.3, min(1.2, lambda_s))
  lambda_t <- max(0.05, min(0.4, lambda_t))
  list(lambda_s=lambda_s, lambda_t=lambda_t, sigma_sp=sigma_sp, sigma_t=sigma_t, sigma_r=sigma_r)
}
