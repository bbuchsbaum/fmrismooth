# minimal neuroim2-friendly array helpers
nv_as_array <- function(x) {
  if (inherits(x, c("NeuroVol","NeuroVec"))) return(as.array(x))
  if (is.array(x)) return(x)
  stop("Expected array or neuroim2::NeuroVol/NeuroVec")
}

nv_wrap_like <- function(proto, arr) {
  if (inherits(proto, "NeuroVec") && requireNamespace("neuroim2", quietly = TRUE)) {
    sp <- tryCatch(neuroim2::space(proto), error = function(e) NULL)
    if (!is.null(sp)) return(neuroim2::NeuroVec(arr, space = sp))
  }
  if (inherits(proto, "NeuroVol") && requireNamespace("neuroim2", quietly = TRUE)) {
    sp <- tryCatch(neuroim2::space(proto), error = function(e) NULL)
    if (!is.null(sp)) return(neuroim2::NeuroVol(arr, space = sp))
  }
  arr
}

assert_numeric_array <- function(x, ndim) {
  if (!is.numeric(x) || length(dim(x)) != ndim) stop("Expected numeric array of dim ", ndim)
  x
}
