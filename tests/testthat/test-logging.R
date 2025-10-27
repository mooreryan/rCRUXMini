describe("logging", {
  it("doesn't get weird in parallel execution", {
    log_file <- tempfile()
    logger <- littlelogger::new_logger(file = log_file)

    with(future::plan(future::multisession, workers = 2), local = TRUE)
    furrr::future_map(1:10, function(x) {
      logger$info("hello %d", x)
    })

    # Read all lines of log_file
    lines <- readLines(log_file)

    # Very basic smoke test to ensure logger lines aren't getting squished
    # together due to parallel file writing.
    expect_equal(length(lines), 10)

    expect_true(all(startsWith(lines, "I, ")))
  })

  it("is magical", {
    log_file <- tempfile()
    logger <- littlelogger::new_logger(file = log_file)

    with(future::plan(future::multisession, workers = 2), local = TRUE)
    furrr::future_map(1:10, function(x) {
      logger$info("hello %d", x)
    })

    # Read all lines of log_file
    lines <- readLines(log_file)

    # Very basic smoke test to ensure logger lines aren't getting squished
    # together due to parallel file writing.
    expect_equal(length(lines), 10)

    expect_true(all(startsWith(lines, "I, ")))
  })
})
