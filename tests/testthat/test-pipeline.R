describe("the pipeline", {
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

    redact_log_details <- function(line) {
      line |>
        stringr::str_replace(
          pattern = "\\d{4}-\\d{2}-\\d{2}",
          replacement = "YYYY-MM-DD"
        ) |>
        stringr::str_replace(
          pattern = "\\d{2}:\\d{2}:\\d{2}\\.\\d{6}",
          replacement = "HH:mm:ss.dddddd"
        ) |>
        stringr::str_replace(
          pattern = "#\\d+",
          replacement = "#PID"
        )
    }

    output_data_files <- output_files |>
      purrr::set_names(output_files_basenames) |>
      purrr::map(function(file) {
        file |>
          readr::read_lines(n_max = 10) |>
          purrr::map(function(line) {
            preview <- line |>
              redact_log_details() |>
              stringr::str_sub(start = 1, end = 70)
            stringr::str_glue(
              "{preview}..."
            )
          })
      })

    expect_snapshot(output_data_files)
    expect_snapshot(result)
  })

  # TODO: this is repeated in multi_primer_multi_db
  it("handles multiple primers and DBs", {
    taxonomy_db <- testthat::test_path(
      "data",
      "small_test_db",
      "taxonomy.db"
    )

    # Run the pipeline with the split DB
    blast_db_multi <- purrr::map(1:10, function(i) {
      testthat::test_path(
        "data",
        "small_test_db",
        sprintf("generated_sequences.part_%03d", i)
      )
    })
    output_directory_multi <- tempfile()
    config_content_multi <- sprintf(
      "forward_primers:
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
reverse_primers:
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
output_directory: '%s'
workers: 2
taxonomy_database: '%s'
blast_databases:
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  - '%s'
  ",
      output_directory_multi,
      taxonomy_db,
      blast_db_multi[[1]],
      blast_db_multi[[2]],
      blast_db_multi[[3]],
      blast_db_multi[[4]],
      blast_db_multi[[5]],
      blast_db_multi[[6]],
      blast_db_multi[[7]],
      blast_db_multi[[8]],
      blast_db_multi[[9]],
      blast_db_multi[[10]]
    )
    config_file_multi <- create_temp_config(config_content_multi)
    result_multi <- pipeline(config = new_config(config_file_multi))

    # Run the pipeline with the same inputs, but on the single DB
    output_directory_single <- tempfile()
    blast_db_single <- testthat::test_path(
      "data",
      "small_test_db",
      "generated_sequences"
    )
    config_content_single <- sprintf(
      "forward_primers:
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
  - 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC'
reverse_primers:
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
  - 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT'
output_directory: '%s'
workers: 2
taxonomy_database: '%s'
blast_databases: '%s'
  ",
      output_directory_single,
      taxonomy_db,
      blast_db_single
    )
    config_file_single <- create_temp_config(config_content_single)
    result_single <- pipeline(config = new_config(config_file_single))

    # Compare the results

    # TODO: check more stuff
    expect_equal(
      result_single$parsed_amplicon_blast_result_distinct_taxonomic_ranks,
      result_multi$parsed_amplicon_blast_result_distinct_taxonomic_ranks
    )

    # Clean up
    unlink(output_directory_multi, recursive = TRUE)
    unlink(config_file_multi)
    unlink(output_directory_single, recursive = TRUE)
    unlink(config_file_single)
  })
})
