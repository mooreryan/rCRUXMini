read_config_file <- function(config_file) {
  yaml::read_yaml(config_file)
}

#' Validate the raw config as read from the YAML file and apply defaults values
#' for any option that the user doesn't specify.
#'
process_config <- function(config) {
  # Validate primer blast parameters if they are provided and set defaults
  config$primer_blast <- new_primer_blast_params(config$primer_blast)

  config$plausible_amplicons <- new_plausible_amplicons_params(
    config$plausible_amplicons
  )
}
