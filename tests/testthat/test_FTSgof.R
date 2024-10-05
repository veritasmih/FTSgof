library(testthat)

test_that("fport_eda works as expected", {
  # Mock data for testing
  set.seed(123)
  mock_data <- matrix(rnorm(5000), nrow = 50, ncol = 100)  # 50 grid points, 20 observations
  
  # Check error handling
  expect_error(fport_eda(mock_data, H = -5), "The parameter 'H' must be a positive integer.")
  expect_error(fport_eda(mock_data, alpha = 1.5), "The 'alpha' parameter must be a value between 0 and 1.")
  
  })

test_that("fport_wn works as expected", {
  # Mock data for testing
  set.seed(123)
  mock_data <- matrix(rnorm(5000), nrow = 50, ncol = 100)  # 50 grid points, 20 observations
  
  # Test if error is thrown for invalid test name
  expect_error(fport_wn(mock_data, test = "invalid_test"), "Please see the documentation for available tests.")
  
})

library(testthat)

test_that("fport_gof handles invalid inputs correctly", {
  f_data <- matrix(rnorm(100), nrow=10, ncol=10)  # Sample functional data
  
  # Test invalid test input
  expect_error(fport_gof(f_data, test = "invalid"), 
               "Please see the documentation for available tests.")
  
  # Test invalid M (non-integer or negative)
  expect_error(fport_gof(f_data, M = -1), 
               "Invalid arguments, M must be a positive integer or NULL.")
  
  # Test invalid f_data (not a matrix)
  expect_error(fport_gof(list(1,2,3), test = "far"), 
               "Invalid arguments, functional data f_data must be passed in matrix form.")
})

test_that("fport_gof handles FAR test correctly", {
  f_data <- matrix(rnorm(5000), nrow=50, ncol=100)  # Sample functional data
  
  # Test FAR(1) model with valid input
  result <- fport_gof(f_data, test = "far", H = 5, M = 10, pplot = FALSE)
  
  # Check that output is a list and contains residuals when residual = TRUE
  result_with_resid <- fport_gof(f_data, test = "far", H = 5, M = 10, residual = TRUE)
  expect_type(result_with_resid, "list")
  expect_true("resid" %in% names(result_with_resid))
})


