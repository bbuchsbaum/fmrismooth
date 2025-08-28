# Test helper functions for fmrismooth package

# Create synthetic 4D fMRI data
create_test_data_4d <- function(dims = c(10, 10, 10, 20), 
                                noise_sd = 1, 
                                signal_mean = 100) {
  signal <- array(signal_mean, dim = dims)
  noise <- array(rnorm(prod(dims), 0, noise_sd), dim = dims)
  signal + noise
}

# Create synthetic 3D volume
create_test_data_3d <- function(dims = c(10, 10, 10), 
                                noise_sd = 1, 
                                signal_mean = 100) {
  signal <- array(signal_mean, dim = dims)
  noise <- array(rnorm(prod(dims), 0, noise_sd), dim = dims)
  signal + noise
}

# Create test mask
create_test_mask <- function(dims, fill = TRUE) {
  if (length(dims) == 4) dims <- dims[1:3]
  array(fill, dim = dims)
}

# Create a spherical mask
create_sphere_mask <- function(dims, center = NULL, radius = NULL) {
  if (length(dims) == 4) dims <- dims[1:3]
  if (is.null(center)) center <- (dims + 1) / 2
  if (is.null(radius)) radius <- min(dims) / 3
  
  mask <- array(FALSE, dim = dims)
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        dist <- sqrt((i - center[1])^2 + (j - center[2])^2 + (k - center[3])^2)
        mask[i, j, k] <- dist <= radius
      }
    }
  }
  mask
}

# Check if arrays are approximately equal
expect_array_equal <- function(x, y, tolerance = 1e-8) {
  expect_equal(dim(x), dim(y))
  expect_true(all(abs(x - y) < tolerance, na.rm = TRUE))
}

# Check dimensions match
expect_dims_equal <- function(x, expected_dims) {
  expect_equal(dim(x), expected_dims)
}

# Check no NaN or Inf values
expect_finite <- function(x) {
  expect_true(all(is.finite(x)))
}