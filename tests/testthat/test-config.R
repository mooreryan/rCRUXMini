describe("creating a new config object from a YAML file", {
  # Setup helper to create temporary config file
  create_temp_config <- function(content) {
    tmp <- tempfile(fileext = ".yaml")
    writeLines(content, tmp)
    tmp
  }

  default_pipeline_config_file <- test_path(
    "data",
    "default_pipeline_config.yml"
  )

  expect_default_params <- function(config) {
    # Defaults for primer blast
    expect_equal(config$primer_blast$evalue, 3e7)
    expect_equal(config$primer_blast$num_alignments, 10000000)
    expect_equal(config$primer_blast$num_threads, 1)
    expect_equal(config$primer_blast$perc_identity, 50)
    expect_equal(config$primer_blast$qcov_hsp_perc, 90)
    expect_equal(config$primer_blast$reward, 2)
    expect_equal(config$primer_blast$task, "blastn-short")
    expect_equal(config$primer_blast$word_size, 7)

    # Defaults for plausible amplicons
    expect_equal(config$plausible_amplicons$minimum_length, 150)
    expect_equal(config$plausible_amplicons$maximum_length, 650)
    expect_equal(config$plausible_amplicons$maximum_mismatches, 4)
    expect_equal(config$plausible_amplicons$ambiguous_run_limit, 5)
  }

  it("reads and parses a valid config file", {
    config <- new_config(default_pipeline_config_file)

    expect_s3_class(config, "rcrux_config")
    expect_s3_class(config$primer_blast, "rcrux_primer_blast_config")
    expect_s3_class(
      config$plausible_amplicons,
      "rcrux_plausible_amplicons_config"
    )
    expect_type(config, "list")
    expect_true("primer_blast" %in% names(config))
    expect_true("plausible_amplicons" %in% names(config))
  })

  it("correctly parses parameters", {
    config <- new_config(default_pipeline_config_file)
    expect_default_params(config)
  })

  it("handles empty or minimal config files with defaults", {
    config_content <- "
primer_blast:

plausible_amplicons:
"

    config_file <- create_temp_config(config_content)
    config <- new_config(config_file)
    expect_default_params(config)

    unlink(config_file)
  })

  it("handles config file with no data at all", {
    config_content <- ""

    config_file <- create_temp_config(config_content)
    config <- new_config(config_file)
    expect_default_params(config)

    unlink(config_file)
  })

  it("handles partial config by using default values for missing params", {
    config_content <- "
primer_blast:
  evalue: 47
  task: 'blastn'

plausible_amplicons:
  minimum_length: 15
  maximum_length: 65
"

    config_file <- create_temp_config(config_content)
    config <- new_config(config_file)

    # Defaults for primer blast
    expect_equal(config$primer_blast$evalue, 47)
    expect_equal(config$primer_blast$num_alignments, 10000000)
    expect_equal(config$primer_blast$num_threads, 1)
    expect_equal(config$primer_blast$perc_identity, 50)
    expect_equal(config$primer_blast$qcov_hsp_perc, 90)
    expect_equal(config$primer_blast$reward, 2)
    expect_equal(config$primer_blast$task, "blastn")
    expect_equal(config$primer_blast$word_size, 7)

    # Defaults for plausible amplicons
    expect_equal(config$plausible_amplicons$minimum_length, 15)
    expect_equal(config$plausible_amplicons$maximum_length, 65)
    expect_equal(config$plausible_amplicons$maximum_mismatches, 4)
    expect_equal(config$plausible_amplicons$ambiguous_run_limit, 5)

    unlink(config_file)
  })

  it("throws an error for non-existent config file", {
    expect_error(
      new_config("non_existent_file.yaml"),
      class = "rcrux_mini_error"
    )
  })
})
