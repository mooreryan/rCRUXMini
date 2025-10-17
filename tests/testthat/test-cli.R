describe("the command line interface", {
  # This is a basic smoke test for the CLI interface for rCRUXMini
  it("works", {
    config_data <- create_config_for_test()
    on.exit(config_data$on_exit(), add = TRUE)

    script_path <- testthat::test_path("../../inst/cli/rCRUXMini.R")
    repo_path <- testthat::test_path("../..")

    run_result <- processx::run(
      "Rscript",
      args = c("--vanilla", script_path, config_data$config_file, repo_path)
    )

    system2("ls", config_data$output_directory, stdout = TRUE) |>
      expect_snapshot()
  })
})
