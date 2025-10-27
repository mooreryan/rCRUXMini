describe("creating a new config object from a YAML file", {
  # Helper to create a minimal valid config with required files/dirs
  create_minimal_valid_config <- function(extra_content = "") {
    # Create temporary files and directories that need to exist
    taxonomy_db <- tempfile(fileext = ".db")
    file.create(taxonomy_db)

    blast_db <- tempfile(fileext = ".nsq")
    file.create(blast_db)
    # Now that we created a fake DB with the nsq suffix, strip it off and use
    # that as the argument.
    blast_db <- stringr::str_replace(blast_db, "\\.nsq$", "")

    output_directory <- tempfile()
    # Don't create output_dir, it must NOT exist

    config_content <- sprintf(
      "forward_primers: 'ACGT'
reverse_primers: 'TGCA'
output_directory: '%s'
taxonomy_database: '%s'
blast_databases: '%s'
%s",
      output_directory,
      taxonomy_db,
      blast_db,
      extra_content
    )

    config_file <- create_temp_config(config_content)

    list(
      config_file = config_file,
      taxonomy_db = taxonomy_db,
      blast_db = blast_db,
      output_directory = output_directory
    )
  }

  expect_default_params <- function(config) {
    # Defaults for general config
    expect_equal(config$query_chunk_count, 1)
    expect_equal(config$ncbi_bin_directory, "")
    expect_equal(config$workers, 1)

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
    temp_files <- create_minimal_valid_config()
    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    config <- new_config(temp_files$config_file)

    expect_type(config, "list")
    expect_s3_class(config, "rcrux_config")
    expect_s3_class(config$primer_blast, "rcrux_primer_blast_config")
    expect_s3_class(
      config$plausible_amplicons,
      "rcrux_plausible_amplicons_config"
    )
  })

  it("handles some non-default top-level params", {
    temp_files <- create_minimal_valid_config(
      "query_chunk_count: 4
ncbi_bin_directory: '/usr/local/ncbi/blast/bin'"
    )
    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    config <- new_config(temp_files$config_file)

    # See the minimal valid config function for some of these values
    expect_equal(config$forward_primers, "ACGT")
    expect_equal(config$reverse_primers, "TGCA")
    expect_equal(config$output_directory, temp_files$output_directory)
    expect_equal(config$taxonomy_database, temp_files$taxonomy_db)
    expect_equal(config$blast_databases, temp_files$blast_db)
    expect_equal(config$query_chunk_count, 4)
    expect_equal(config$ncbi_bin_directory, "/usr/local/ncbi/blast/bin")
  })

  it("correctly parses multiple primers and blast databases", {
    temp_files <- create_minimal_valid_config()
    blast_db2 <- tempfile(fileext = ".nsq")
    file.create(blast_db2)
    blast_db2 <- stringr::str_replace(blast_db2, "\\.nsq$", "")

    config_content <- sprintf(
      "forward_primers:
  - 'ACGT'
  - 'AAAA'
reverse_primers:
  - 'TGCA'
  - 'TTTT'
output_directory: '%s'
taxonomy_database: '%s'
blast_databases:
  - '%s'
  - '%s'",
      temp_files$output_dir,
      temp_files$taxonomy_db,
      temp_files$blast_db,
      blast_db2
    )

    config_file <- create_temp_config(config_content)

    on.exit({
      unlink(config_file)
      unlink(temp_files$taxonomy_db)
      unlink(stringr::str_glue("{temp_files$blast_db}.nsq"))
      unlink(stringr::str_glue(
        "{blast_db2}.nsq"
      ))
    })

    config <- new_config(config_file)

    expect_equal(config$forward_primers, c("ACGT", "AAAA"))
    expect_equal(config$reverse_primers, c("TGCA", "TTTT"))
    expect_equal(config$blast_databases, c(temp_files$blast_db, blast_db2))
  })

  it("applies default values for optional parameters", {
    temp_files <- create_minimal_valid_config()
    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    config <- new_config(temp_files$config_file)
    expect_default_params(config)
  })

  it("handles partial config by using default values for missing params", {
    temp_files <- create_minimal_valid_config(
      "query_chunk_count: 8
workers: 3

primer_blast:
  evalue: 47
  task: 'blastn'

plausible_amplicons:
  minimum_length: 15
  maximum_length: 65"
    )

    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    config <- new_config(temp_files$config_file)

    # General params (custom and defaults)
    expect_equal(config$query_chunk_count, 8)
    expect_equal(config$workers, 3)
    expect_equal(config$ncbi_bin_directory, "")

    # Primer blast params (custom and defaults)
    expect_equal(config$primer_blast$evalue, 47)
    expect_equal(config$primer_blast$num_alignments, 10000000)
    expect_equal(config$primer_blast$num_threads, 1)
    expect_equal(config$primer_blast$perc_identity, 50)
    expect_equal(config$primer_blast$qcov_hsp_perc, 90)
    expect_equal(config$primer_blast$reward, 2)
    expect_equal(config$primer_blast$task, "blastn")
    expect_equal(config$primer_blast$word_size, 7)

    # Plausible amplicons params (custom and defaults)
    expect_equal(config$plausible_amplicons$minimum_length, 15)
    expect_equal(config$plausible_amplicons$maximum_length, 65)
    expect_equal(config$plausible_amplicons$maximum_mismatches, 4)
    expect_equal(config$plausible_amplicons$ambiguous_run_limit, 5)
  })

  it("validates that the config file exists", {
    expect_error(
      new_config("non_existent_file.yaml"),
      class = "rcrux_mini_error",
      regexp = "does not exist"
    )
  })

  it("validates that output_directory does not exist", {
    temp_files <- create_minimal_valid_config()
    # Create the output directory (which should fail validation)
    dir.create(temp_files$output_dir)

    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
      unlink(temp_files$output_dir, recursive = TRUE)
    })

    expect_error(
      new_config(temp_files$config_file),
      class = "rcrux_mini_error",
      regexp = "already exists"
    )
  })

  it("validates that taxonomy_db_path exists", {
    temp_files <- create_minimal_valid_config()
    # TODO: use the add param for all on.exit calls
    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    # Delete the taxonomy DB file, which should trigger the validation error
    unlink(temp_files$taxonomy_db)

    expect_error(
      new_config(temp_files$config_file),
      class = "rcrux_mini_error",
      regexp = "does not exist"
    )
  })

  it("validates that blast_db_paths exist", {
    temp_files <- create_minimal_valid_config()

    # Delete the taxonomy DB file, which should trigger the validation error
    unlink(stringr::str_glue("{temp_files$blast_db}.nsq"))

    # TODO: use the add param for all on.exit calls
    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(stringr::str_glue("{temp_files$blast_db}.nsq"))
    })

    expect_error(
      new_config(temp_files$config_file),
      # TODO: it would be better for this to be an rcrux error or even a
      # SnailBLAST error
      class = "checkmateError",
      regexp = "No .nsq or .psq file found"
    )
  })

  it("validates that ncbi_bin_directory exists when provided", {
    # We don't create this file!
    ncbi_bin <- tempfile()

    temp_files <- create_minimal_valid_config(
      sprintf("ncbi_bin_directory: '%s'", ncbi_bin)
    )

    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
    })

    expect_error(
      new_config(temp_files$config_file),
      class = "rcrux_mini_error",
      regexp = "does not exist"
    )
  })

  it("accepts ncbi_bin_directory when it exists", {
    ncbi_bin <- tempfile()
    dir.create(ncbi_bin)

    temp_files <- create_minimal_valid_config(
      sprintf("ncbi_bin_directory: '%s'", ncbi_bin)
    )

    on.exit({
      unlink(temp_files$config_file)
      unlink(temp_files$taxonomy_db)
      unlink(temp_files$blast_db)
      unlink(ncbi_bin, recursive = TRUE)
    })

    config <- new_config(temp_files$config_file)
    expect_equal(config$ncbi_bin_directory, ncbi_bin)
  })

  it("throws an error for non-existent config file", {
    expect_error(
      new_config("non_existent_file.yaml"),
      class = "rcrux_mini_error"
    )
  })

  it("throws an error when required general parameters are missing", {
    config_content <- ""

    config_file <- create_temp_config(config_content)

    on.exit(unlink(config_file))

    expect_error(
      new_config(config_file),
      # TODO: would be nice if this were an rcrux error
      class = "checkmateError",
      regexp = "Names must include"
    )
  })
})
