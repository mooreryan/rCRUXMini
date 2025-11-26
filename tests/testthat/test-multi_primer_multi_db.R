describe("Multiple primers and multiple DBs", {
  # 10 blast DBs
  # 3 workers
  # 7 query chunk counts
  it("it works with non-divisible numbers", {
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
query_chunk_count: 7
workers: 3
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
