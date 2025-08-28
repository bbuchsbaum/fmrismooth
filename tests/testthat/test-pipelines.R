test_that("fmrismooth_default with auto parameters", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 15), noise_sd = 3, signal_mean = 100)
  
  # Run with auto parameters and no robust stage
  result <- fmrismooth_default(
    vec,
    robust = "none",
    final = "joint_bilateral",
    auto_params = TRUE,
    tr = 2.0,
    target_fwhm_mm = 5
  )
  
  # Check dimensions preserved
  expect_equal(dim(result), dim(vec))
  
  # Check finite
  expect_finite(result)
})

test_that("fmrismooth_default with robust smoothing", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 15), noise_sd = 3, signal_mean = 100)
  
  # Test with Huber robust
  result_huber <- fmrismooth_default(
    vec,
    robust = "huber",
    final = "joint_bilateral",
    auto_params = TRUE
  )
  
  expect_equal(dim(result_huber), dim(vec))
  expect_finite(result_huber)
  
  # Test with Tukey robust
  result_tukey <- fmrismooth_default(
    vec,
    robust = "tukey",
    final = "joint_bilateral",
    auto_params = TRUE
  )
  
  expect_equal(dim(result_tukey), dim(vec))
  expect_finite(result_tukey)
})

test_that("fmrismooth_default with guided filter final stage", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 15), noise_sd = 3)
  
  # Use guided filter (implemented as weak bilateral)
  result <- fmrismooth_default(
    vec,
    robust = "none",
    final = "guided_filter",
    auto_params = TRUE
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fmrismooth_default with manual parameters", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 15))
  
  # Use manual parameters (auto_params = FALSE)
  result <- fmrismooth_default(
    vec,
    robust = "none",
    final = "joint_bilateral",
    auto_params = FALSE  # Will use default manual params
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fmrismooth_default with T1 guide", {
  set.seed(123)
  
  # Create test data and T1 guide
  dims <- c(8, 8, 8, 15)
  vec <- create_test_data_4d(dims = dims)
  t1 <- create_test_data_3d(dims = dims[1:3], signal_mean = 150)
  
  # Run with T1 guide
  result <- fmrismooth_default(
    vec,
    robust = "none",
    t1 = t1,
    auto_params = TRUE
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fmrismooth_mppca_pipeline works", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 15), noise_sd = 3)
  
  # Run MP-PCA pipeline
  result <- fmrismooth_mppca_pipeline(
    vec,
    sigma_mode = "global",
    sigma_sp = 2.5,
    sigma_t = 0.5,
    sigma_r = 15,
    lattice_blur_iters = 1L
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fmrismooth_mppca_pipeline with T1 and tissue probs", {
  set.seed(123)
  
  # Create test data
  dims <- c(8, 8, 8, 15)
  vec <- create_test_data_4d(dims = dims)
  t1 <- create_test_data_3d(dims = dims[1:3], signal_mean = 150)
  
  # Create tissue probability maps (simplified)
  prob1 <- create_test_data_3d(dims = dims[1:3], signal_mean = 0.5)
  prob1 <- pmax(0, pmin(1, prob1 / max(prob1)))  # Normalize to [0,1]
  
  # Run with guides
  result <- fmrismooth_mppca_pipeline(
    vec,
    t1 = t1,
    probs = list(prob1),
    sigma_mode = "patch"
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("pipelines reduce noise variance", {
  set.seed(123)
  
  # Create noisy data
  dims <- c(8, 8, 8, 15)
  clean <- array(100, dim = dims)
  noise <- array(rnorm(prod(dims), 0, 5), dim = dims)
  noisy <- clean + noise
  
  # Test default pipeline
  result1 <- fmrismooth_default(noisy, robust = "none", auto_params = TRUE)
  var_original <- var(as.vector(noisy))
  var_smoothed1 <- var(as.vector(result1))
  expect_true(var_smoothed1 < var_original)
  
  # Test MP-PCA pipeline
  result2 <- fmrismooth_mppca_pipeline(noisy, sigma_mode = "global")
  var_smoothed2 <- var(as.vector(result2))
  expect_true(var_smoothed2 < var_original)
})