# Setup helper to create temporary config file
create_temp_config <- function(content) {
  tmp <- tempfile(fileext = ".yaml")
  writeLines(content, tmp)
  tmp
}

# Helper to create a minimal valid config with required files/dirs
create_config_for_test <- function(
  forward_primers = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC",
  reverse_primers = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT"
) {
  taxonomy_db <- testthat::test_path(
    "data",
    "small_test_db",
    "taxonomy.db"
  )

  blast_db <- testthat::test_path(
    "data",
    "small_test_db",
    "generated_sequences"
  )

  output_directory <- tempfile()
  # We don't create output_dir because it must NOT exist

  config_content <- sprintf(
    "forward_primers: '%s'
reverse_primers: '%s'
output_directory: '%s'
taxonomy_database: '%s'
blast_databases: '%s'
workers: 2
",
    forward_primers,
    reverse_primers,
    output_directory,
    taxonomy_db,
    blast_db
  )

  # TODO: this function is from the other test
  config_file <- create_temp_config(config_content)

  on_exit <- function() {
    unlink(output_directory, recursive = TRUE)
    unlink(config_file)
  }

  list(
    config_file = config_file,
    taxonomy_db = taxonomy_db,
    blast_db = blast_db,
    output_directory = output_directory,
    on_exit = on_exit
  )
}


redact_tmp_filename <- function(filename, replacement = "REDACTED") {
  # Replace hexadecimal strings (temp file identifiers)
  # Pattern matches: underscore + hex chars + file extension
  gsub("_[0-9a-f]+\\.", paste0("_", replacement, "."), filename)
}
