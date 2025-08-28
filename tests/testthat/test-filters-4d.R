test_that("fast_bilateral_lattice4d preserves dimensions", {
  set.seed(123)
  
  # Create 4D test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10), noise_sd = 2)
  
  # Apply bilateral filter
  result <- fast_bilateral_lattice4d(
    vec, 
    sigma_sp = 2.0, 
    sigma_t = 0.5, 
    sigma_r = 10.0,
    blur_iters = 1L
  )
  
  # Check dimensions preserved
  expect_equal(dim(result), dim(vec))
  
  # Check finite values
  expect_finite(result)
})

test_that("fast_bilateral_lattice4d with spatial guide", {
  set.seed(123)
  
  # Create data and guide
  dims <- c(8, 8, 8, 10)
  vec <- create_test_data_4d(dims = dims)
  guide_spatial <- create_test_data_3d(dims = dims[1:3], signal_mean = 50)
  
  # Apply with guide
  result <- fast_bilateral_lattice4d(
    vec,
    sigma_sp = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0,
    guide_spatial = guide_spatial
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fast_bilateral_lattice4d with design matrix", {
  set.seed(123)
  
  # Create data
  dims <- c(8, 8, 8, 10)
  vec <- create_test_data_4d(dims = dims)
  
  # Create simple design vector
  design <- seq(0, 1, length.out = dims[4])
  
  # Apply with design
  result <- fast_bilateral_lattice4d(
    vec,
    sigma_sp = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0,
    design = design,
    sigma_d = 1.0
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("fast_bilateral_joint4d works correctly", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10), noise_sd = 3)
  guide <- apply(vec, 1:3, mean)  # Mean across time as guide
  
  # Apply joint bilateral
  result <- fast_bilateral_joint4d(
    vec, 
    guide_spatial = guide,
    sigma_s = 2.0,
    sigma_t = 0.5,
    sigma_r = 10.0,
    passes = 1L
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("st_guided_filter4d basic functionality", {
  set.seed(123)
  
  # Create 4D test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10), noise_sd = 3)
  
  # Apply guided filter
  result <- st_guided_filter4d(
    vec,
    radius_s = 2L,
    radius_t = 1L,
    eps_s = 0.5^2,
    eps_t = 0.3^2
  )
  
  # Check dimensions
  expect_equal(dim(result), dim(vec))
  
  # Check finite
  expect_finite(result)
})

test_that("mp_pca4d denoises correctly", {
  set.seed(123)
  
  # Create noisy 4D data
  dims <- c(10, 10, 10, 20)
  clean <- array(100, dim = dims)
  noise <- array(rnorm(prod(dims), 0, 3), dim = dims)
  noisy <- clean + noise
  
  # Apply MP-PCA
  denoised <- mp_pca4d(
    noisy,
    sigma_mode = "global",
    patch = c(5L, 5L, 5L),
    tw = 10L,
    stride = c(3L, 3L, 3L, 5L)
  )
  
  # Check dimensions
  expect_equal(dim(denoised), dim(noisy))
  
  # Check finite
  expect_finite(denoised)
  
  # Variance should be reduced
  var_original <- var(as.vector(noisy))
  var_denoised <- var(as.vector(denoised))
  expect_true(var_denoised < var_original)
})

test_that("stv_denoise4d works", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10), noise_sd = 2)
  
  # Apply STV denoising
  result <- stv_denoise4d(
    vec,
    lambda_s = 0.5,
    lambda_t = 0.2,
    iters = 10L
  )
  
  expect_equal(dim(result), dim(vec))
  expect_finite(result)
})

test_that("stv_robust4d with different loss functions", {
  set.seed(123)
  
  # Create test data
  vec <- create_test_data_4d(dims = c(8, 8, 8, 10), noise_sd = 2)
  
  # Test Huber loss
  result_huber <- stv_robust4d(
    vec,
    loss = "huber",
    lambda_s = 0.5,
    lambda_t = 0.2
  )
  
  expect_equal(dim(result_huber), dim(vec))
  expect_finite(result_huber)
  
  # Test Tukey loss
  result_tukey <- stv_robust4d(
    vec,
    loss = "tukey",
    lambda_s = 0.5,
    lambda_t = 0.2
  )
  
  expect_equal(dim(result_tukey), dim(vec))
  expect_finite(result_tukey)
})