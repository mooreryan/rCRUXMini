# This is essentially the same as test-pipeline, except that it writes the fully
# serialized value. The reason it is not in the regular test, is that the
# serialized output makes the snapshot diffing too tedious. So do _not_ change
# this test very often!

test_that("pipeline works", {
  config_data <- create_config_for_test()
  on.exit(config_data$on_exit(), add = TRUE)

  config <- new_config(config_data$config_file)
  result <- pipeline(config = config)

  # And this is the full data structure
  expect_snapshot_value(result, style = "serialize")
})
