test_that("VST forward and inverse are inverses above noise floor", {
  # Test data (exclude values below sqrt(2)*sigma where forward clamps to 0)
  sigma <- 5
  x <- seq(10, 100, by = 10)
  
  z <- vst_forward_rician(x, sigma)
  x_recovered <- vst_inverse_rician(z, sigma)
  expect_equal(x_recovered, x, tolerance = 1e-10)
})

test_that("VST forward handles edge cases", {
  sigma <- 5
  
  # Zero input
  expect_equal(vst_forward_rician(0, sigma), 0)
  
  # Negative values (should be clipped to 0 by pmax)
  result <- vst_forward_rician(-10, sigma)
  expect_equal(result, 0)
  
  # Large values
  x <- 1000
  result <- vst_forward_rician(x, sigma)
  expect_true(result > 0)
  expect_true(is.finite(result))
})

test_that("VST inverse handles edge cases", {
  sigma <- 5
  
  # Zero input
  result <- vst_inverse_rician(0, sigma)
  expect_equal(result, sqrt(2 * sigma^2))
  
  # Negative values (should still work)
  result <- vst_inverse_rician(-10, sigma)
  expect_equal(result, sqrt(2 * sigma^2))
  
  # Large values
  z <- 1000
  result <- vst_inverse_rician(z, sigma)
  expect_true(result > 0)
  expect_true(is.finite(result))
})

test_that("vst_wrap applies denoising correctly", {
  set.seed(123)
  
  # Create test data
  data_4d <- create_test_data_4d(dims = c(5, 5, 5, 10), noise_sd = 2, signal_mean = 100)
  
  # Simple identity denoiser for testing
  identity_denoise <- function(x, ...) x
  
  # Apply VST wrap
  result <- vst_wrap(data_4d, sigma = 2, denoise_fun = identity_denoise)
  
  # Should preserve dimensions
  expect_equal(dim(result), dim(data_4d))
  
  # Should be finite
  expect_finite(result)
  
  # Test with actual smoothing
  smooth_denoise <- function(x, ...) {
    # Simple averaging along time
    apply(x, 1:3, mean)
  }
  
  # This reduces to 3D; wrapper should error because shape must be preserved
  expect_error(vst_wrap(data_4d, sigma = 2, denoise_fun = smooth_denoise))
})

test_that("vst_wrap estimates sigma when not provided", {
  set.seed(123)
  
  # Create 4D test data
  data_4d <- create_test_data_4d(dims = c(5, 5, 5, 10), noise_sd = 3, signal_mean = 100)
  
  # Identity denoiser
  identity_denoise <- function(x, ...) x
  
  # Should estimate sigma automatically for 4D data
  result <- vst_wrap(data_4d, sigma = NULL, denoise_fun = identity_denoise)
  expect_equal(dim(result), dim(data_4d))
  expect_finite(result)
  
  # Should error for 3D data without sigma
  data_3d <- create_test_data_3d(dims = c(5, 5, 5))
  expect_error(vst_wrap(data_3d, sigma = NULL, denoise_fun = identity_denoise), 
               "sigma must be supplied for 3D")
})
