# TODO: rather than use the example pipeline config, write a config file specific for the test (START HERE)

describe("the pipeline", {
  # TODO: this is duplicated from test-config

  it("works", {
    config_data <- create_config_for_test()
    on.exit(config_data$on_exit(), add = TRUE)

    config <- new_config(config_data$config_file)
    result <- pipeline(config = config)

    system2("ls", config_data$output_directory, stdout = TRUE) |>
      expect_snapshot()

    output_files <- stringr::str_glue("{config_data$output_directory}/*") |>
      Sys.glob() |>
      purrr::discard(function(file) {
        stringr::str_ends(file, pattern = ".rds")
      })

    output_files_basenames <- output_files |> purrr::map(basename)

    output_data_files <- output_files |>
      purrr::set_names(output_files_basenames) |>
      purrr::map(function(file) {
        file |>
          readr::read_lines(n_max = 10) |>
          purrr::map(function(line) {
            preview <- line |> stringr::str_sub(start = 1, end = 70)
            stringr::str_glue(
              "{preview}..."
            )
          })
      })

    expect_snapshot(output_data_files)
    expect_snapshot(result)
  })
})
