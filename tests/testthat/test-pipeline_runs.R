# This is essentially testing nothing. It is here in cases where you just need
# to run the pipeline and either manually examine the results using `browser()`
# or to see if it crashes.

test_that("the pipeline runs", {
  config_data <- create_config_for_test()
  on.exit(config_data$on_exit(), add = TRUE)

  config <- new_config(config_data$config_file)
  result <- pipeline(config = config)

  testthat::expect_true(!is.null(result))
})
