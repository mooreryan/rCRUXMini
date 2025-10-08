new_config <- function(config_file) {
  config_file |>
    .read_config_file() |>
    .validate_config_and_apply_default_params() |>
    structure(class = "rcrux_pipeline_config")
}

.read_config_file <- function(config_file) {
  if (!file.exists(config_file)) {
    abort_rcrux_mini_error("Config file does not exist")
  }

  yaml::read_yaml(config_file)
}

#' Validate the raw config as read from the YAML file and apply defaults values
#' for any option that the user doesn't specify.
#'
.validate_config_and_apply_default_params <- function(config) {
  # Validate primer blast parameters if they are provided and set defaults
  config$primer_blast <- new_primer_blast_params(config$primer_blast)

  config$plausible_amplicons <- new_plausible_amplicons_params(
    config$plausible_amplicons
  )

  config
}

# Basic parsing functions are needed because users frequently pass blast params
# as scientific notation, and not all the forms are recognized by the yaml
# parser. E.g., 3e7 will parse to a character rather than a numeric.
#
# Ensure that you call these on all params that should parse to either integers
# or numerics in their respective parse/validation functions.

.parse_integer <- function(x) {
  as.integer(x)
}

.parse_numeric <- function(x) {
  as.numeric(x)
}
