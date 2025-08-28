test_that("recommend_params returns sensible defaults", {
  set.seed(123)
  
  # Create test data
  data_4d <- create_test_data_4d(dims = c(10, 10, 10, 20), noise_sd = 2, signal_mean = 100)
  
  # Get recommendations
  params <- recommend_params(data_4d, mask = NULL, tr = 2.0, target_fwhm_mm = 5)
  
  # Check that all expected parameters are returned
  expect_named(params, c("lambda_s", "lambda_t", "sigma_sp", "sigma_t", "sigma_r"))
  
  # Check parameters are within reasonable bounds
  expect_true(params$lambda_s >= 0.3 && params$lambda_s <= 1.2)
  expect_true(params$lambda_t >= 0.05 && params$lambda_t <= 0.4)
  expect_true(params$sigma_sp >= 1.0)
  expect_true(params$sigma_t >= 0.3 && params$sigma_t <= 0.8)
  expect_true(params$sigma_r >= 5)
  
  # All should be finite
  expect_true(all(is.finite(unlist(params))))
})

test_that("recommend_params handles different TR values", {
  set.seed(123)
  data_4d <- create_test_data_4d(dims = c(10, 10, 10, 20))
  
  # Short TR
  params_short <- recommend_params(data_4d, tr = 0.5, target_fwhm_mm = 5)
  
  # Long TR
  params_long <- recommend_params(data_4d, tr = 4.0, target_fwhm_mm = 5)
  
  # Temporal sigma should adapt to TR
  expect_true(params_short$sigma_t != params_long$sigma_t)
  
  # Lambda_t should also adapt
  expect_true(params_short$lambda_t != params_long$lambda_t)
})

test_that("recommend_params handles different target smoothness", {
  set.seed(123)
  data_4d <- create_test_data_4d(dims = c(10, 10, 10, 20))
  
  # Small smoothing
  params_small <- recommend_params(data_4d, target_fwhm_mm = 2)
  
  # Large smoothing
  params_large <- recommend_params(data_4d, target_fwhm_mm = 10)
  
  # Spatial sigma should scale with target FWHM
  expect_true(params_small$sigma_sp < params_large$sigma_sp)
})

test_that("recommend_params works with mask", {
  set.seed(123)
  dims <- c(10, 10, 10, 20)
  data_4d <- create_test_data_4d(dims = dims)
  mask <- create_sphere_mask(dims)
  
  # Should work with mask
  params <- recommend_params(data_4d, mask = mask)
  
  expect_named(params, c("lambda_s", "lambda_t", "sigma_sp", "sigma_t", "sigma_r"))
  expect_true(all(is.finite(unlist(params))))
})

test_that("recommend_params accepts pre-computed sigma", {
  set.seed(123)
  data_4d <- create_test_data_4d(dims = c(10, 10, 10, 20))
  
  # With pre-computed sigma
  params <- recommend_params(data_4d, sigma_mppca = 2.5)
  
  # Should still work
  expect_named(params, c("lambda_s", "lambda_t", "sigma_sp", "sigma_t", "sigma_r"))
  expect_true(all(is.finite(unlist(params))))
})