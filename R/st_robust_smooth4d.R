#' Spatiotemporal robust smoothing (Huber data term + anisotropic TV)
#'
#' Robust fMRI denoising that downweights outlier frames/voxels (Huber M-estimation)
#' while preserving edges via anisotropic space-time TV. Solved by outer IRLS loops
#' (to update Huber data weights) with an inner Chambolle-Pock primal-dual solver.
#'
#' @param vec 4D fMRI (array or neuroim2::NeuroVec)
#' @param lambda_s spatial TV weight
#' @param lambda_t temporal TV weight
#' @param delta Huber threshold (in signal units) for the data term; smaller => more robust
#' @param outer_loops number of outer IRLS updates (2-5 typical)
#' @param iters number of inner primal-dual iterations per outer loop (20-40 typical)
#' @param mask 3D mask (outside voxels passed through unchanged)
#' @param motion_params optional T x 6 (or 3) realignment params to reduce temporal smoothing during high motion
#' @param sigma_m motion weighting scale (same units as framewise displacement proxy)
#' @return 4D array or NeuroVec
#' @export
st_robust_smooth4d <- function(vec,
                               lambda_s = 0.8,
                               lambda_t = 0.25,
                               delta    = 1.5,
                               outer_loops = 3L,
                               iters       = 25L,
                               mask = NULL,
                               motion_params = NULL,
                               sigma_m = 0.5) {
  X <- nv_as_array(vec); X <- assert_numeric_array(X, 4L)
  d <- dim(X)
  # align / validate mask
  m <- tryCatch(align_mask_3d(mask, vec), error=function(e) nv_check_mask(mask, d))

  # temporal weights from motion (optional)
  wt <- rep(1.0, d[4])
  if (!is.null(motion_params)) {
    mp <- as.matrix(motion_params)
    if (!ncol(mp) %in% c(3,6)) stop("motion_params should have 6 (or 3) cols.")
    if (nrow(mp) != d[4]) stop("motion_params rows must equal number of frames.")
    dmp <- sqrt(rowSums(rbind(rep(0, ncol(mp)), diff(mp))^2))
    wt  <- as.numeric(exp(-(dmp^2) / (2 * sigma_m^2)))
  }

  Y <- tv_robust_4d_cpp(X, m,
                        lambda_s, lambda_t,
                        delta,
                        outer_loops, iters,
                        wt)
  nv_wrap_like(vec, Y)
}
