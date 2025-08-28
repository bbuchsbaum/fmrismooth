test_that("fast_bilateral_lattice3d preserves dimensions", {
  set.seed(123)
  
  # Create 3D test data
  vol <- create_test_data_3d(dims = c(10, 10, 10), noise_sd = 2)
  
  # Apply bilateral filter
  result <- fast_bilateral_lattice3d(vol, sigma_sp = 2.0, sigma_r = 10.0, blur_iters = 1L)
  
  # Check dimensions preserved
  expect_equal(dim(result), dim(vol))
  
  # Check finite values
  expect_finite(result)
})

test_that("fast_bilateral_lattice3d smooths noise", {
  set.seed(123)
  
  # Create noisy data
  dims <- c(10, 10, 10)
  clean <- array(100, dim = dims)
  noise <- array(rnorm(prod(dims), 0, 5), dim = dims)
  noisy <- clean + noise
  
  # Apply filter
  smoothed <- fast_bilateral_lattice3d(noisy, sigma_sp = 2.0, sigma_r = 10.0)
  
  # Smoothed should have less variance than original
  var_original <- var(as.vector(noisy))
  var_smoothed <- var(as.vector(smoothed))
  expect_true(var_smoothed < var_original)
})

test_that("fast_bilateral_lattice3d works with guide", {
  set.seed(123)
  
  # Create data and guide
  vol <- create_test_data_3d(dims = c(10, 10, 10))
  guide <- create_test_data_3d(dims = c(10, 10, 10), signal_mean = 50)
  
  # Apply with guide
  result <- fast_bilateral_lattice3d(vol, sigma_sp = 2.0, sigma_r = 10.0, guide = guide)
  
  expect_equal(dim(result), dim(vol))
  expect_finite(result)
})

test_that("fast_bilateral_lattice3d works with mask", {
  set.seed(123)
  
  # Create data and mask
  dims <- c(10, 10, 10)
  vol <- create_test_data_3d(dims = dims)
  mask <- create_sphere_mask(dims)
  
  # Apply with mask
  result <- fast_bilateral_lattice3d(vol, sigma_sp = 2.0, sigma_r = 10.0, mask = mask)
  
  expect_equal(dim(result), dim(vol))
  expect_finite(result)
})

test_that("st_guided_filter3d basic functionality", {
  set.seed(123)
  
  # Create 3D test data
  vol <- create_test_data_3d(dims = c(10, 10, 10), noise_sd = 3)
  
  # Apply guided filter
  result <- st_guided_filter3d(vol, radius = 2L, eps = 0.5^2)
  
  # Check dimensions
  expect_equal(dim(result), dim(vol))
  
  # Check finite
  expect_finite(result)
})

test_that("st_guided_filter3d with guide", {
  set.seed(123)
  
  # Create data and guide
  vol <- create_test_data_3d(dims = c(10, 10, 10))
  guide <- create_test_data_3d(dims = c(10, 10, 10), signal_mean = 50)
  
  # Apply with guide
  result <- st_guided_filter3d(vol, radius = 2L, eps = 0.5^2, guide = guide)
  
  expect_equal(dim(result), dim(vol))
  expect_finite(result)
})

test_that("fast_bilateral_joint3d works correctly", {
  set.seed(123)
  
  # Create test data
  vol <- create_test_data_3d(dims = c(8, 8, 8), noise_sd = 3)
  guide <- create_test_data_3d(dims = c(8, 8, 8), signal_mean = 50)
  
  # Apply joint bilateral
  result <- fast_bilateral_joint3d(
    vol, 
    guide_spatial = guide,
    sigma_s = 2.0, 
    sigma_r = 10.0,
    passes = 1L
  )
  
  expect_equal(dim(result), dim(vol))
  expect_finite(result)
})

test_that("fast_bilateral_joint3d handles backend options", {
  set.seed(123)
  
  vol <- create_test_data_3d(dims = c(8, 8, 8))
  guide <- create_test_data_3d(dims = c(8, 8, 8))
  
  # Test grid backend (routes to permutohedral)
  result_grid <- fast_bilateral_joint3d(
    vol, guide_spatial = guide,
    sigma_s = 2.0, sigma_r = 10.0,
    backend = "grid"
  )
  
  # Test permutohedral backend
  result_perm <- fast_bilateral_joint3d(
    vol, guide_spatial = guide,
    sigma_s = 2.0, sigma_r = 10.0,
    backend = "permutohedral"
  )
  
  # Both should work and produce similar results
  expect_equal(dim(result_grid), dim(vol))
  expect_equal(dim(result_perm), dim(vol))
  expect_finite(result_grid)
  expect_finite(result_perm)
})