#' Create an empty tibble from a readr column specification
#'
#' @param col_spec A readr column specification object (cols object)
#' @return A tibble with 0 rows and columns matching the specification
#'
#' @examples
#' # Create a column specification
#' spec <- cols(
#'   id = col_integer(),
#'   name = col_character(),
#'   score = col_double(),
#'   date = col_date(),
#'   active = col_logical()
#' )
#'
#' # Generate empty tibble
#' empty_df <- empty_tibble_from_col_spec(spec)
#' print(empty_df)
#'
empty_tibble_from_col_spec <- function(col_spec) {
  checkmate::assert_class(col_spec, "col_spec")
  checkmate::assert_named(col_spec$cols)

  # Create empty vectors for each column based on its type
  columns <- purrr::map(col_spec$cols, function(col_type) {
    class_name <- class(col_type)[1]

    if (class_name == "collector_guess") {
      abort_rcrux_mini_error("you aren't allowed to use collector_guess")
    }

    switch(
      class_name,
      "collector_logical" = logical(0),
      "collector_integer" = integer(0),
      "collector_double" = numeric(0),
      "collector_number" = numeric(0),
      "collector_character" = character(0),
      "collector_datetime" = as.POSIXct(character(0)),
      "collector_date" = as.Date(character(0)),
      "collector_time" = hms::as_hms(character(0)),
      "collector_factor" = {
        if (!is.null(col_type$levels)) {
          factor(character(0), levels = col_type$levels)
        } else {
          factor(character(0))
        }
      },
      # Skip columns are not included
      "collector_skip" = NULL,
      # Default to character for unknown types
      character(0)
    )
  })

  # Remove NULL entries (skipped columns)
  non_empty_columns <- columns[!purrr::map_lgl(columns, is.null)]
  non_empty_col_names <- names(col_spec$cols)[!purrr::map_lgl(columns, is.null)]

  # Create tibble
  result <- tibble::tibble(!!!non_empty_columns)
  colnames(result) <- non_empty_col_names
  result
}
