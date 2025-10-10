describe("the pipeline", {
  # TODO: create a config file on the fly
  config_file <- test_path(
    "data",
    "example_pipeline_config.yml"
  )

  it("works", {
    tmpdir_top <- tempdir()

    output_directory_path <- file.path(tmpdir_top, "rcrux_output")
    on.exit(unlink(output_directory_path, recursive = TRUE), add = TRUE)

    config <- new_config(config_file)

    result <- pipeline(
      forward_primers = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC",
      reverse_primers = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT",
      output_directory_path = output_directory_path,
      blast_db_paths = test_path(
        "data",
        "small_test_db",
        "generated_sequences"
      ),
      taxonomy_db_path = test_path(
        "data",
        "small_test_db",
        "taxonomy.db"
      ),
      query_chunk_count = 1,
      ncbi_bin_directory = NULL,
      config = config
    )

    system2("ls", output_directory_path, stdout = TRUE) |>
      expect_snapshot()

    output_files <- stringr::str_glue("{output_directory_path}/*") |>
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
