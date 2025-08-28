#' @keywords internal
nv_as_array <- function(x) {
  if (is.array(x)) {
    # Ensure double storage explicitly
    y <- array(as.numeric(x), dim = dim(x))
    storage.mode(y) <- "double"
    return(y)
  }
  cls <- class(x)
  if (any(cls %in% c("NeuroVol", "NeuroVec"))) {
    aa <- tryCatch(as.array(x), error = function(e) NULL)
    if (is.array(aa)) { y <- array(as.numeric(aa), dim = dim(aa)); storage.mode(y) <- "double"; return(y) }
    aa <- tryCatch(x[], error = function(e) NULL)
    if (is.array(aa)) { y <- array(as.numeric(aa), dim = dim(aa)); storage.mode(y) <- "double"; return(y) }
  }
  stop("Expected a numeric array or a neuroim2 NeuroVol/NeuroVec.")
}

#' @keywords internal
nv_wrap_like <- function(template, arr) {
  stopifnot(is.array(arr))
  cls <- class(template)
  if (requireNamespace("neuroim2", quietly = TRUE) &&
      any(cls %in% c("NeuroVol","NeuroVec"))) {
    sp <- tryCatch(neuroim2::space(template), error = function(e) NULL)
    if (!is.null(sp)) {
      if (length(dim(arr)) == 3L && inherits(template, "NeuroVol")) {
        return(neuroim2::NeuroVol(arr, space = sp))
      }
      if (length(dim(arr)) == 4L && inherits(template, "NeuroVec")) {
        return(neuroim2::NeuroVec(arr, space = sp))
      }
    }
  }
  arr
}

#' @keywords internal
nv_check_mask <- function(mask, dims) {
  if (is.null(mask)) {
    if (length(dims) == 3L) return(array(TRUE, dims))
    if (length(dims) == 4L) return(array(TRUE, dims[1:3]))
  }
  # allow NeuroVol or array mask
  m <- if (any(class(mask) == "NeuroVol")) tryCatch(as.array(mask), error=function(e) NULL) else mask
  if (is.null(m)) m <- nv_as_array(mask)
  if (length(dim(m)) == 4L && identical(dim(m)[1:3], dims[1:3])) {
    if (dim(m)[4] == dims[length(dims)]) {
      m <- apply(m, 1:3, function(v) any(v != 0))
    }
  }
  # Accept integer vs double differences in dims
  td <- if (length(dims) == 4L) dims[1:3] else dims
  if (!(length(dim(m)) == 3L && all(as.integer(dim(m)) == as.integer(td)))) {
    stop("mask grid != target grid; expected ", paste(td, collapse = "x"))
  }
  m != 0
}

#' @keywords internal
assert_numeric_array <- function(x, nd) {
  if (!is.array(x) || length(dim(x)) != nd) stop("Expected a numeric ", nd, "D array.")
  y <- array(as.numeric(x), dim = dim(x)); storage.mode(y) <- "double"
  y
}
