#' 4D space-time TV (ROF) denoising via Chambolle-Pock
#'
#' Performs anisotropic total variation denoising in space and time.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param lambda_s Spatial TV weight.
#' @param lambda_t Temporal TV weight.
#' @param iters Number of primal-dual iterations (positive integer).
#' @param mask Optional 3D logical/0-1 mask.
#' @param tau Algorithm primal step size; if `NULL`, a sensible default is used.
#' @param sigma Algorithm dual step size; if `NULL`, a sensible default is used.
#' @param theta Over-relaxation parameter.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
stv_denoise4d <- function(vec,
                          lambda_s = 0.8, lambda_t = 0.2,
                          iters = 30L, mask = NULL,
                          tau = NULL, sigma = NULL, theta = 1.0) {
  if (!is.numeric(lambda_s) || length(lambda_s)!=1L || !is.finite(lambda_s)) stop("lambda_s must be a finite numeric scalar")
  if (!is.numeric(lambda_t) || length(lambda_t)!=1L || !is.finite(lambda_t)) stop("lambda_t must be a finite numeric scalar")
  if (!is.numeric(iters) || length(iters)!=1L || is.na(iters) || as.integer(iters) <= 0L) stop("iters must be a positive integer")
  X <- nv_as_array(vec); X <- assert_numeric_array(X, 4L)
  d <- dim(X); m <- nv_check_mask(mask, d)
  if (is.null(tau) || is.null(sigma)) {
    tau   <- 1.0 / sqrt(8.0)
    sigma <- 1.0 / sqrt(8.0)
  }
  Y <- .Call("_fmrismooth_tv_denoise_4d_cpp", X, m,
             as.numeric(lambda_s), as.numeric(lambda_t),
             as.integer(iters), as.numeric(tau),
             as.numeric(sigma), as.numeric(theta),
             PACKAGE = "fmrismooth")
  nv_wrap_like(vec, Y)
}
