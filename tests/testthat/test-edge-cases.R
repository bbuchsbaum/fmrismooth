test_that("functions handle NULL mask correctly", {
  set.seed(123)
  
  # 3D data with NULL mask
  vol <- create_test_data_3d(dims = c(8, 8, 8))
  result_3d <- fast_bilateral_lattice3d(vol, sigma_sp = 2.0, sigma_r = 10.0, mask = NULL)
  expect_equal(dim(result_3d), dim(vol))
  
  # 4D data with NULL mask
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10))
  result_4d <- fast_bilateral_lattice4d(vec, sigma_sp = 2.0, sigma_t = 0.5, sigma_r = 10.0, mask = NULL)
  expect_equal(dim(result_4d), dim(vec))
})

test_that("functions error on wrong dimensions", {
  # 2D data to 3D function
  data_2d <- array(1:100, dim = c(10, 10))
  expect_error(fast_bilateral_lattice3d(data_2d))
  
  # 3D data to 4D function
  data_3d <- array(1:1000, dim = c(10, 10, 10))
  expect_error(fast_bilateral_lattice4d(data_3d))
  
  # 5D data
  data_5d <- array(1:10000, dim = c(10, 10, 10, 10, 1))
  expect_error(mp_pca4d(data_5d))
})

test_that("functions handle single time frame gracefully", {
  # Single frame should error for temporal operations
  single_frame <- create_test_data_4d(dims = c(8, 8, 8, 1))
  
  # estimate_sigma_rician needs at least 3 frames
  expect_error(estimate_sigma_rician(single_frame), "at least 3 frames")
})

test_that("functions handle extreme parameter values", {
  set.seed(123)
  vec <- create_test_data_4d(dims = c(5, 5, 5, 10))
  
  # Very large sigma (should still work)
  result <- fast_bilateral_lattice4d(
    vec,
    sigma_sp = 100.0,  # Very large
    sigma_t = 10.0,    # Very large
    sigma_r = 1000.0   # Very large
  )
  expect_finite(result)
  
  # Very small sigma (should still work)
  result <- fast_bilateral_lattice4d(
    vec,
    sigma_sp = 0.1,   # Very small
    sigma_t = 0.01,   # Very small
    sigma_r = 0.1     # Very small
  )
  expect_finite(result)
})

test_that("functions handle mismatched guide dimensions", {
  # 4D data with wrong guide dimensions
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10))
  wrong_guide <- create_test_data_3d(dims = c(5, 5, 5))  # Wrong size
  
  # Should error when guides don't match
  expect_error(fast_bilateral_joint4d(
    vec,
    guide_spatial = wrong_guide,
    sigma_s = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0
  ))
})

test_that("functions handle empty/all-false masks", {
  set.seed(123)
  
  # Create data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10))
  
  # All-false mask
  empty_mask <- array(FALSE, dim = c(8, 8, 8))
  
  # Should still run but may produce warnings or specific behavior
  result <- fast_bilateral_lattice4d(
    vec,
    sigma_sp = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0,
    mask = empty_mask
  )
  
  expect_equal(dim(result), dim(vec))
})

test_that("functions handle NaN and Inf in input", {
  # Data with NaN
  vec_nan <- create_test_data_4d(dims = c(5, 5, 5, 10))
  vec_nan[1, 1, 1, 1] <- NaN
  
  # Most functions should handle NaN gracefully
  result <- fast_bilateral_lattice4d(
    vec_nan,
    sigma_sp = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0
  )
  expect_equal(dim(result), dim(vec_nan))
  
  # Data with Inf
  vec_inf <- create_test_data_4d(dims = c(5, 5, 5, 10))
  vec_inf[1, 1, 1, 1] <- Inf
  
  result <- fast_bilateral_lattice4d(
    vec_inf,
    sigma_sp = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0
  )
  expect_equal(dim(result), dim(vec_inf))
})

test_that("pipelines handle missing optional dependencies gracefully", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10))
  
  # Should work even without neuroim2 package
  result <- fmrismooth_default(
    vec,
    robust = "none",
    auto_params = TRUE
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("functions validate parameter types", {
  vec <- create_test_data_4d(dims = c(5, 5, 5, 10))
  
  # Non-numeric sigma
  expect_error(fast_bilateral_lattice4d(
    vec,
    sigma_sp = "two",  # Should be numeric
    sigma_t = 0.5,
    sigma_r = 10.0
  ))
  
  # Negative iterations
  expect_error(stv_denoise4d(
    vec,
    lambda_s = 0.5,
    lambda_t = 0.2,
    iters = -5L  # Should be positive
  ))
})