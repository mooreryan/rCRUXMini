cli_args <- commandArgs(trailingOnly = TRUE)

num_args <- length(cli_args)

if (!(num_args %in% c(1, 2))) {
  message("\n\n")
  message(
    "usage: Rscript --vanilla pipeline.R </path/to/config.yml> [path/to/rCRUXMini]"
  )
  message(
    "  - The 1st argument is the path to a YAML file containing the parameters."
  )
  message(
    "  - The 2nd argument is an optional path to a local rCRUXMini repository."
  )
  message(
    "    Use it if you want to run a local (development) version of rCRUX."
  )
  quit(save = "no", status = 1)
}

if (num_args == 2) {
  # TODO: switch this out to not use devtools so we don't have to make users
  # install it
  devtools::load_all(cli_args[[2]], quiet = TRUE)
} else {
  library(rCRUXMini, quietly = TRUE)
}

# TODO: add to imports
config <- yaml::read_yaml(cli_args[[1]])

print(config)

# Capture the result here.  This prevents prenting NULL to the console when
# calling Rscript.
result <- do.call(what = run, args = config)
