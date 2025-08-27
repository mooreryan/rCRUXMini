test_that("empty_tibble_from_col_spec handles all readr col types", {
  result <- empty_tibble_from_col_spec(readr::cols(
    lgl = readr::col_logical(),
    int = readr::col_integer(),
    dbl = readr::col_double(),
    num = readr::col_number(),
    chr = readr::col_character(),
    dttm = readr::col_datetime(),
    date = readr::col_date(),
    time = readr::col_time(),
    fct = readr::col_factor(levels = c("A", "B")),
    skip = readr::col_skip()
  ))

  testthat::expect_equal(
    object = result,
    expected = tibble::tibble(
      lgl = logical(0),
      int = integer(0),
      dbl = double(0),
      num = double(0),
      chr = character(0),
      dttm = as.POSIXct(character(0)),
      date = as.Date(character(0)),
      time = hms::hms(numeric(0)),
      fct = factor(character(0), levels = c("A", "B"))
      # note: no "skip" column expected
    )
  )
})

test_that("empty_tibble_from_col_spec validates input type", {
  # Not a col_spec
  input <- list(cols = list())

  testthat::expect_error(
    empty_tibble_from_col_spec(input),
    class = "checkmateError"
  )
})

test_that("empty_tibble_from_col_spec requires named columns in col_spec", {
  testthat::expect_error(
    empty_tibble_from_col_spec(readr::cols(readr::col_integer())),
    class = "checkmateError"
  )
})

test_that("empty_tibble_from_col_spec errors on collector_guess", {
  guess_spec <- readr::cols(x = readr::col_guess())
  testthat::expect_error(
    empty_tibble_from_col_spec(guess_spec),
    class = "rcrux_mini_error"
  )
})

test_that("empty_tibble_from_col_spec preserves order and drops skipped columns", {
  spec <- readr::cols(
    a = readr::col_integer(),
    drop_me = readr::col_skip(),
    b = readr::col_character()
  )
  result <- empty_tibble_from_col_spec(spec)

  checkmate::expect_data_frame(
    result,
    nrow = 0,
    types = c("integer", "character")
  )
  checkmate::expect_names(names(result), identical.to = c("a", "b"))
})

test_that("empty_tibble_from_col_spec returns 0-column tibble when all columns skipped", {
  spec <- readr::cols(only_skip = readr::col_skip())
  result <- empty_tibble_from_col_spec(spec)

  checkmate::expect_data_frame(result, nrows = 0, ncols = 0)
})

test_that("empty_tibble_from_col_spec handles factor without levels", {
  spec <- readr::cols(fct = readr::col_factor())
  result <- empty_tibble_from_col_spec(spec)

  testthat::expect_true(is.factor(result$fct))
  testthat::expect_equal(length(levels(result$fct)), 0)
})

test_that("empty_tibble_from_col_spec falls back to character for unknown collectors", {
  # Manually craft a col_spec with an unknown collector class
  unknown_collector <- structure(list(), class = "collector_foo")

  spec <- structure(
    list(cols = list(unknown = unknown_collector)),
    class = "col_spec"
  )

  result <- empty_tibble_from_col_spec(spec)

  checkmate::expect_data_frame(result, nrows = 0, types = "character")
  checkmate::expect_names(colnames(result), identical.to = "unknown")
})
