#' One-liner fMRI smoother with auto-params, robust stage, and final-stage choice
#'
#' Convenience pipeline that optionally applies a robust space-time smoothing
#' stage, followed by a joint bilateral or guided-filter final stage. Parameters
#' can be estimated automatically from data and voxel size.
#'
#' @param vec 4D fMRI array or `NeuroVec`.
#' @param robust Robust loss for pre-smoothing: one of `"huber"`, `"tukey"`, or `"none"`.
#' @param final Final stage method: `"joint_bilateral"` or `"guided_filter"`.
#' @param backend Joint bilateral backend label (currently routes to permutohedral).
#' @param t1 Optional 3D anatomical guide.
#' @param probs Optional list of 3D probability maps used as additional guides.
#' @param mask Optional 3D logical/0-1 mask.
#' @param interp_guide Guide interpolation mode when resampling with neuroim2.
#' @param interp_mask Mask interpolation mode when resampling with neuroim2.
#' @param auto_params If `TRUE`, auto-estimate parameters via `recommend_params`.
#' @param tr Repetition time (seconds) for auto-params; optional.
#' @param target_fwhm_mm Target spatial FWHM (in mm) for auto-params.
#' @param motion_params Optional motion regressors for temporal weighting.
#' @param design Optional design regressor for design-aware joint bilateral.
#' @return Smoothed 4D data, wrapped like `vec` when possible.
#' @export
fmrismooth_default <- function(
  vec,
  robust = c("huber","tukey","none"),
  final  = c("joint_bilateral","guided_filter"),
  backend = c("grid","permutohedral"),
  t1 = NULL, probs = NULL, mask = NULL,
  interp_guide = 1L, interp_mask = 0L,
  auto_params = TRUE, tr = NULL, target_fwhm_mm = 5,
  # robust weights:
  motion_params = NULL, design = NULL
) {
  robust <- match.arg(robust); final <- match.arg(final); backend <- match.arg(backend)
  X <- nv_as_array(vec); d <- dim(X)
  mask <- tryCatch(align_mask_3d(mask, vec, interp = interp_mask), error=function(e) array(TRUE, d[1:3]))
  rec <- if (auto_params) recommend_params(vec, mask = mask, tr = tr, target_fwhm_mm = target_fwhm_mm)
          else list(lambda_s=0.8, lambda_t=0.2, sigma_sp=2.5, sigma_t=0.5, sigma_r=15)

  Y <- vec
  if (robust != "none") {
    Y <- stv_robust4d(
      Y, loss = robust,
      lambda_s = rec$lambda_s, lambda_t = rec$lambda_t,
      motion_params = motion_params, sigma_m = 0.5,
      mask = mask
    )
  }
  guide_spatial <- if (is.null(t1)) apply(nv_as_array(Y), 1:3, mean) else t1
  if (final == "joint_bilateral") {
    Y <- fast_bilateral_joint4d(
      Y, guide_spatial = guide_spatial,
      sigma_s = rec$sigma_sp, sigma_t = rec$sigma_t, sigma_r = rec$sigma_r,
      passes = 1L, mask = mask,
      interp_guide = interp_guide, interp_mask = interp_mask,
      backend = backend, design = design, guides = probs
    )
  } else {
    # minimal guided filter fallback: run joint bilateral with large sigma_r (weak guidance)
    Y <- fast_bilateral_joint4d(
      Y, guide_spatial = guide_spatial,
      sigma_s = rec$sigma_sp, sigma_t = rec$sigma_t, sigma_r = max(50, rec$sigma_r*4),
      passes = 1L, mask = mask,
      interp_guide = interp_guide, interp_mask = interp_mask,
      backend = backend, design = design, guides = probs
    )
  }
  nv_wrap_like(vec, Y)
}
