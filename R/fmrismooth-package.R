#' fmrismooth: Fast Edge-Preserving 3D/4D Smoothing and Denoising for fMRI
#'
#' Edge-preserving guided filters, MP-PCA, TV denoising, and lattice bilateral filters,
#' tailored to fMRI volumes/time-series. If 'neuroim2' is available, guides/masks are
#' automatically resampled to the fMRI grid.
#'
#' @useDynLib fmrismooth, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
