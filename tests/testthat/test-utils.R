test_that("nv_as_array converts arrays correctly", {
  # Test with plain array
  arr <- array(1:8, dim = c(2, 2, 2))
  result <- fmrismooth:::nv_as_array(arr)
  expect_equal(result, arr)
  expect_true(is.numeric(result))
  
  # Test with non-array input (should error)
  expect_error(fmrismooth:::nv_as_array(list(1, 2, 3)))
  expect_error(fmrismooth:::nv_as_array(data.frame(a = 1:3)))
})

test_that("nv_wrap_like preserves array when no neuroim2", {
  # Without neuroim2, should return array unchanged
  template <- array(0, dim = c(10, 10, 10))
  arr <- array(1:1000, dim = c(10, 10, 10))
  result <- fmrismooth:::nv_wrap_like(template, arr)
  expect_equal(result, arr)
})

test_that("nv_check_mask handles various inputs", {
  # NULL mask should create all-TRUE mask
  dims <- c(10, 10, 10, 20)
  mask <- fmrismooth:::nv_check_mask(NULL, dims)
  expect_equal(dim(mask), dims[1:3])
  expect_true(all(mask))
  
  # 3D mask for 4D data
  mask_3d <- array(TRUE, dim = dims[1:3])
  result <- fmrismooth:::nv_check_mask(mask_3d, dims)
  expect_equal(dim(result), dims[1:3])
  
  # Logical conversion
  mask_numeric <- array(1, dim = dims[1:3])
  result <- fmrismooth:::nv_check_mask(mask_numeric, dims)
  expect_type(result, "logical")
})

test_that("assert_numeric_array validates correctly", {
  # Valid array
  arr <- array(1:27, dim = c(3, 3, 3))
  result <- fmrismooth:::assert_numeric_array(arr, 3)
  expect_true(is.numeric(result))
  
  # Wrong dimensions
  expect_error(fmrismooth:::assert_numeric_array(arr, 2))
  expect_error(fmrismooth:::assert_numeric_array(arr, 4))
  
  # Non-array
  expect_error(fmrismooth:::assert_numeric_array(1:10, 1))
})

test_that("estimate_sigma_rician estimates noise correctly", {
  set.seed(123)
  # Create 4D data with known noise
  true_sigma <- 2.0
  dims <- c(10, 10, 10, 20)
  clean_data <- array(100, dim = dims)
  
  # Add Rician noise (simplified - just Gaussian for testing)
  noise <- array(rnorm(prod(dims), 0, true_sigma), dim = dims)
  noisy_data <- sqrt((clean_data + noise)^2)
  
  # Estimate sigma
  estimated <- estimate_sigma_rician(noisy_data)
  
  # Should be positive
  expect_true(estimated > 0)
  
  # Should be finite
  expect_true(is.finite(estimated))
  
  # Should be roughly in the right ballpark (within 50% of true value)
  expect_true(estimated > true_sigma * 0.5)
  expect_true(estimated < true_sigma * 1.5)
})
