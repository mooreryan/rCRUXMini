test_that("pipeline works", {
  tmpdir_top <- tempdir()

  output_directory_path <- file.path(tmpdir_top, "rcrux_output")
  on.exit(unlink(output_directory_path, recursive = TRUE), add = TRUE)

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
    # TODO: this needs to change to the run length (see find_amplicons)
    wildcards = "NNNNNNNNNNNN"
  )

  # This is for a human readable preview
  expect_snapshot(result)
  # And this is the full data structure
  expect_snapshot_value(result, style = "serialize")
})
