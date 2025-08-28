#' Robust space-time smoothing (Huber or Tukey) + anisotropic TV
#'
#' Applies a robust data term (Huber or Tukey) with spatial/temporal TV priors,
#' optionally with motion-aware temporal weights.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param lambda_s Spatial TV weight.
#' @param lambda_t Temporal TV weight.
#' @param loss Robust loss: `'huber'` or `'tukey'`.
#' @param delta Huber threshold; if `NULL`, derived from noise.
#' @param cthresh Tukey threshold; if `NULL`, derived from noise.
#' @param alpha Tukey shape parameter.
#' @param temporal_weights Optional length-(T or T-1) temporal edge weights.
#' @param motion_params Optional motion regressors to derive temporal weights.
#' @param sigma_m Motion weight scale.
#' @param iters Iterations for the robust solver.
#' @param mask Optional 3D logical/0-1 mask.
#' @param tau Algorithm primal step size; if `NULL`, a sensible default is used.
#' @param sigma Algorithm dual step size; if `NULL`, a sensible default is used.
#' @param theta Over-relaxation parameter.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
stv_robust4d <- function(vec,
                         lambda_s = 0.8, lambda_t = 0.2,
                         loss = c('huber','tukey'),
                         delta = NULL, cthresh = NULL,
                         alpha = 1.0,
                         temporal_weights = NULL,
                         motion_params = NULL, sigma_m = 0.5,
                         iters = 35L, mask = NULL,
                         tau = NULL, sigma = NULL, theta = 1.0) {
  X <- nv_as_array(vec); X <- assert_numeric_array(X, 4L)
  d <- dim(X)
  m <- tryCatch(align_mask_3d(mask, vec), error=function(e) array(TRUE, d[1:3]))
  loss <- match.arg(loss)
  if (is.null(tau) || is.null(sigma)) { tau <- 1.0/sqrt(8.0); sigma <- 1.0/sqrt(8.0) }
  if (loss=='huber' && is.null(delta)) {
    sig <- estimate_sigma_rician(X, mask = m); delta <- 1.345 * sig
  }
  if (loss=='tukey' && is.null(cthresh)) {
    sig <- estimate_sigma_rician(X, mask = m); cthresh <- 4.685 * sig
  }
  wt_edge <- NULL
  if (!is.null(temporal_weights)) {
    tw <- as.numeric(temporal_weights)
    if (length(tw) == d[4]) wt_edge <- pmin(tw[-d[4]], tw[-1]) else if (length(tw) == d[4]-1) wt_edge <- tw
    wt_edge[!is.finite(wt_edge)] <- 1; wt_edge <- pmax(0, pmin(1, wt_edge))
  } else if (!is.null(motion_params)) {
    mp <- as.matrix(motion_params)
    if (nrow(mp) != d[4]) stop("motion_params rows must equal T")
    dmp <- sqrt(rowSums(rbind(rep(0, ncol(mp)), diff(mp))^2))
    wf <- exp(-(dmp^2) / (2 * sigma_m^2))
    wt_edge <- pmin(wf[-d[4]], wf[-1])
  }
  if (loss=='huber') {
    # Note: tv_robust_4d_cpp expects different parameters (delta, outer_loops, iters, wt_in)
    # We need to adapt the call
    outer_loops <- as.integer(3)  # default outer loops
    wt_in <- rep(1.0, d[4])  # default temporal weights
    if (!is.null(wt_edge)) wt_in[-1] <- wt_edge
    
    Y <- tv_robust_4d_cpp(X, m,
                          lambda_s, lambda_t, 
                          delta, outer_loops, as.integer(iters),
                          wt_in)
  } else {
    Y <- tv_tukey_4d_cpp(X, m,
                         lambda_s, lambda_t,
                         cthresh, alpha,
                         if (is.null(wt_edge)) numeric(0) else as.numeric(wt_edge),
                         as.integer(iters), tau, sigma, theta)
  }
  nv_wrap_like(vec, Y)
}
