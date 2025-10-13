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

  it("handles multiple primers and DBs", {
    taxonomy_db <- testthat::test_path(
      "data",
      "small_test_db",
      "taxonomy.db"
    )

    blast_db_1 <- testthat::test_path(
      "data",
      "small_test_db",
      "generated_sequences.part_001"
    )

    blast_db_2 <- testthat::test_path(
      "data",
      "small_test_db",
      "generated_sequences.part_002"
    )

    blast_db_3 <- testthat::test_path(
      "data",
      "small_test_db",
      "generated_sequences.part_003"
    )

    output_directory <- tempfile()
    # We don't create output_dir because it must NOT exist

    config_content <- sprintf(
      "forward_primers:
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
  reverse_primers:
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
  output_directory: '%s'
  taxonomy_database: '%s'
  blast_databases:
  - '%s'
  - '%s'
  - '%s'
  ",
      output_directory,
      taxonomy_db,
      blast_db_1,
      blast_db_2,
      blast_db_3
    )

    # TODO: this function is from the other test
    config_file <- create_temp_config(config_content)

    on_exit <- function() {
      unlink(output_directory, recursive = TRUE)
      unlink(config_file)
    }

    config_data <- create_config_for_test()
    on.exit(config_data$on_exit(), add = TRUE)

    config <- new_config(config_data$config_file)
    result <- pipeline(config = config)
  })
})
