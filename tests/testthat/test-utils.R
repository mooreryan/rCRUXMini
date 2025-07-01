library(testthat)
library(checkmate)

test_that("enumerate_ambiguity handles standard DNA bases", {
  result <- enumerate_ambiguity("ATCG")
  expect_equal(result, "ATCG")
  expect_length(result, 1)
})

test_that("enumerate_ambiguity expands single ambiguity codes correctly", {
  # M -> A or C
  result <- enumerate_ambiguity("M")
  expect_equal(sort(result), c("A", "C"))
  expect_length(result, 2)

  # N represents any of the four bases
  result <- enumerate_ambiguity("N")
  expect_equal(sort(result), c("A", "C", "G", "T"))
  expect_length(result, 4)
})

test_that("enumerate_ambiguity expands multiple ambiguities correctly", {
  # M -> A C
  # W -> A T
  result <- enumerate_ambiguity("MW")
  expected <- c("AA", "CA", "AT", "CT")
  expect_equal(sort(result), sort(expected))
  expect_length(result, 4)
})

test_that("enumerate_ambiguity handles mixed standard and ambiguous bases", {
  # M -> A C
  result <- enumerate_ambiguity("AMT")
  expected <- c("AAT", "ACT")
  expect_equal(sort(result), sort(expected))
  expect_length(result, 2)
})

test_that("enumerate_ambiguity generates correct number of combinations", {
  # Each ambiguity code multiplies the possibilities
  # MR = 2 * 2 = 4 combinations
  result <- enumerate_ambiguity("MR")
  expect_length(result, 4)

  # MRW = 2 * 2 * 2 = 8 combinations
  result <- enumerate_ambiguity("MRW")
  expect_length(result, 8)

  # KV = 2 * 3 = 6 combinations
  result <- enumerate_ambiguity("KV")
  expect_length(result, 6)

  # BN = 3 * 4 = 12 combinations
  result <- enumerate_ambiguity("BN")
  expect_length(result, 12)
})

test_that("enumerate_ambiguity handles edge cases", {
  # Empty string
  result <- enumerate_ambiguity("")
  expect_equal(result, "")
  expect_length(result, 1)

  # Unknown characters should be passed through
  result <- enumerate_ambiguity("X")
  expect_equal(result, "X")
  expect_length(result, 1)

  # Validates input
  expect_error(enumerate_ambiguity(letters), "length")
  expect_error(enumerate_ambiguity(234), "type")
})

test_that("enumerate_ambiguity produces valid DNA sequences", {
  result <- enumerate_ambiguity("MRWS")

  # All results should be same length as input
  expect_true(all(nchar(result) == 4))

  # All results should only contain A, T, C, G
  valid_chars <- "^[ATCG]+$"
  grepl(pattern = valid_chars, x = result) |>
    all() |>
    expect_true()

  # Should have 2*2*2*2 = 16 combinations
  expect_length(result, 16)

  # All combinations should be unique
  expect_equal(length(result), length(unique(result)))
})

test_that("enumerate_ambiguities", {
  result <- enumerate_ambiguities("AN")
  expected <- c("AA", "AC", "AG", "AT")

  expect_equal(sort(result), expected)
})

test_that("enumerate_ambiguities", {
  result <- enumerate_ambiguities(c("AN", "NA"))
  expected <- sort(c(
    # From AN
    "AA",
    "AC",
    "AG",
    "AT",
    # From NA
    "AA",
    "CA",
    "GA",
    "TA"
  ))

  expect_equal(sort(result), expected)
})


test_that("run_command validates inputs correctly", {
  skip_if(Sys.info()["sysname"] == "Windows")

  # Test invalid command type
  expect_error(run_command(123), "type")

  expect_error(run_command(c("echo", "test")), "length")

  # Test invalid args type
  expect_error(run_command("echo", 123), "type")

  expect_error(run_command("echo", list("test")), "type")
})

test_that("run_command executes simple commands successfully", {
  skip_if(Sys.info()["sysname"] == "Windows")

  result <- run_command("echo", "Hello")
  expected_output <- "Hello"

  expect_list(result)
  expect_names(
    names(result),
    must.include = c(
      "exit_code",
      "stdout",
      "stderr",
      "warnings",
      "full_command"
    )
  )

  expect_equal(result$exit_code, 0)
  expect_equal(trimws(result$stdout), expected_output)
  expect_equal(result$stderr, "")
  expect_equal(result$full_command, "echo Hello")
})

test_that("run_command handles NULL arguments", {
  skip_if(Sys.info()["sysname"] == "Windows")

  result <- run_command("echo", args = NULL)

  expect_equal(result$exit_code, 0)
  expect_equal(result$stdout, "")
  expect_equal(result$stderr, "")
  expect_equal(result$full_command, "echo")
})

test_that("run_command handles multiple arguments", {
  skip_if(Sys.info()["sysname"] == "Windows")

  result <- run_command("echo", c("Hello", "World"))

  expect_equal(result$exit_code, 0)
  expect_equal(trimws(result$stdout), "Hello World")
  expect_equal(result$stderr, "")
  expect_equal(result$full_command, "echo Hello World")
})

test_that("run_command handles commands with zero exit codes", {
  skip_if(Sys.info()["sysname"] == "Windows")

  # Arrange
  script_file <- tempfile()
  on.exit({
    unlink(script_file)
  })
  script <- 'printf "hello"; printf "world" >&2; exit 0'
  write(x = script, file = script_file)

  # Act
  result <- run_command("sh", script_file)

  # Assert
  expect_equal(result$exit_code, 0)
  expect_equal(result$stdout, "hello")
  expect_equal(result$stderr, "world")
})


test_that("run_command handles commands with non-zero exit codes", {
  skip_if(Sys.info()["sysname"] == "Windows")

  # Arrange
  script_file <- tempfile()
  on.exit({
    unlink(script_file)
  })
  script <- 'printf "hello"; printf "world" >&2; exit 1'
  write(x = script, file = script_file)

  # Act
  result <- run_command("sh", script_file)

  # Assert
  expect_equal(result$exit_code, 1)
  expect_equal(result$stdout, "hello")
  expect_equal(result$stderr, "world")
})

test_that("run_command handles empty output correctly", {
  skip_if(Sys.info()["sysname"] == "Windows")

  # Unix command that does nothing
  result <- run_command("true")

  expect_equal(result$exit_code, 0)
  expect_equal(result$stdout, "")
  expect_equal(result$stderr, "")
})

test_that("run_command respects timeout parameter", {
  skip_if(Sys.info()["sysname"] == "Windows")

  # This should timeout
  start_time <- Sys.time()
  result <- run_command("sleep", "5", timeout = 1)
  end_time <- Sys.time()

  # The command should have been terminated due to timeout
  expect_true(result$exit_code != 0)
  expect_true(as.numeric(end_time - start_time) < 2)
  expect_true(stringr::str_length(result$warnings) > 0)
})

test_that("run_command handles commands that don't exist", {
  skip_if(Sys.info()["sysname"] == "Windows")

  result <- run_command("nonexistent_command_12345")

  # Should return non-zero exit code
  expect_true(result$exit_code != 0)
  expect_equal(result$stdout, "")
  expect_true(stringr::str_length(result$warnings) > 0)
})

test_that("run_command return structure is consistent", {
  skip_if(Sys.info()["sysname"] == "Windows")

  result <- run_command("echo", "test")

  expect_type(result, "list")
  expect_true(length(result) >= 5)
  expect_names(
    names(result),
    must.include = c(
      "exit_code",
      "stdout",
      "stderr",
      "warnings",
      "full_command"
    )
  )

  expect_int(result$exit_code)
  expect_string(result$stdout)
  expect_string(result$stderr)
  expect_character(result$warnings)
  expect_string(result$full_command)
})
