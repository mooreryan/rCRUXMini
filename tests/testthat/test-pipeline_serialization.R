# This is essentially the same as test-pipeline, except that it writes the fully
# serialized value. The reason it is not in the regular test, is that the
# serialized output makes the snapshot diffing too tedious. So do _not_ change
# this test very often!

test_that("pipeline works", {
  tmpdir_top <- tempdir()

  output_directory_path <- file.path(tmpdir_top, "rcrux_output")
  on.exit(unlink(output_directory_path, recursive = TRUE), add = TRUE)

  config <- new_config(test_path(
    "data",
    "default_pipeline_config.yml"
  ))

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

  # And this is the full data structure
  expect_snapshot_value(result, style = "serialize")
})
