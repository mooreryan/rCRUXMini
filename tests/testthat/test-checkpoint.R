describe("checkpointing", {
  write_empty_tmpfile <- function() {
    file <- tempfile()
    cat(NULL, file = file)
    file
  }

  it("runs the expression if the checkfile doesn't exist", {
    file <- tempfile()

    result <- checkpoint(file, function() {
      "hello, world!"
    })

    expect_equal(result, "hello, world!")
  })

  it("does NOT run the expr if the checkfile already exists", {
    file <- write_empty_tmpfile()

    result <- checkpoint(file, function() {
      "hello, world!"
    })
    testthat::expect_null(result)
  })
})
