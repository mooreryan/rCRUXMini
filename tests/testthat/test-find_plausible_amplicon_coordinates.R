# The little ascii art diagrams aren't to scale, treat them as relative
# positions.

describe("calculating product length for primers in correct orientation", {
  it("calculates length for normal forward orientation (F before R)", {
    #  F ---------->                        <------------ R  Primers
    # ====================================================== DNA Target
    #    ---------->************************<------------    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 100,
      forward_stop = 120,
      reverse_start = 500,
      reverse_stop = 480
    )

    expect_equal(product_length, 400)
  })

  it("calculates length for reverse orientation (R before F)", {
    #  R ---------->                        <------------ F  Primers
    # ====================================================== DNA Target
    #    ---------->************************<------------    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 500,
      forward_stop = 480,
      reverse_start = 100,
      reverse_stop = 120
    )

    expect_equal(product_length, 400)
  })

  it("handles abutting primers (F before R)", {
    #  F ----------><------------ R  Primers
    # ============================== DNA Target
    #    ----------><------------    Product
    #    ^         ^
    #    100       120
    #               ^           ^
    #               121         140
    product_length <- .calculate_product_length(
      forward_start = 100,
      forward_stop = 120,
      reverse_start = 140,
      reverse_stop = 121
    )

    expect_equal(product_length, 40)
  })

  it("handles abutting primers (R before F)", {
    #  R ----------><------------ F  Primers
    # ============================== DNA Target
    #    ----------><------------    Product
    #    ^         ^
    #    100       120
    #               ^           ^
    #               121         140
    product_length <- .calculate_product_length(
      forward_start = 140,
      forward_stop = 121,
      reverse_start = 100,
      reverse_stop = 120
    )

    expect_equal(product_length, 40)
  })

  it("overlapping primers in correct orientation are okay (F before R)", {
    #  F ---------->
    #         <------------ R
    # ======================== DNA Target
    #    ^         ^
    #    100       150
    #         ^            ^
    #         125          175
    product_length <- .calculate_product_length(
      forward_start = 100,
      forward_stop = 150,
      reverse_start = 175,
      reverse_stop = 125
    )

    expect_equal(product_length, 75)
  })

  it("overlapping primers in correct orientation are okay (R before F)", {
    #  R ---------->
    #         <------------ F
    # ======================== DNA Target
    #    ^         ^
    #    100       150
    #         ^            ^
    #         125          175
    product_length <- .calculate_product_length(
      forward_start = 175,
      forward_stop = 125,
      reverse_start = 100,
      reverse_stop = 150
    )

    expect_equal(product_length, 75)
  })

  it("completly coverlapping primers in correct orientation are okay (F before R)", {
    #  F ---------->
    #    <---------- R
    # ======================== DNA Target
    #    ^         ^
    #    100       150
    product_length <- .calculate_product_length(
      forward_start = 100,
      forward_stop = 150,
      reverse_start = 150,
      reverse_stop = 100
    )

    expect_equal(product_length, 50)
  })

  it("completly coverlapping primers in correct orientation are okay (R before F)", {
    #  R ---------->
    #    <---------- F
    # ======================== DNA Target
    #    ^         ^
    #    100       150
    product_length <- .calculate_product_length(
      forward_start = 150,
      forward_stop = 100,
      reverse_start = 100,
      reverse_stop = 150
    )

    expect_equal(product_length, 50)
  })

  it("works with vectorized inputs", {
    product_lengths <- .calculate_product_length(
      forward_start = c(100, 200, 300),
      forward_stop = c(120, 230, 340),
      reverse_start = c(200, 301, 402),
      reverse_stop = c(180, 270, 360)
    )

    expect_equal(product_lengths, c(100, 101, 102))
  })
})

describe("calculating product length for primers in invalid orientations", {
  it("returns NA with both primers like this: ---> (F before R)", {
    #  F ---------->                        ------------> R  Primers
    # ====================================================== DNA Target
    #    ---------->************************------------>    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 100,
      forward_stop = 120,
      reverse_start = 480,
      reverse_stop = 500
    )

    expect_true(is.na(product_length))
  })

  it("returns NA with both primers like this: ---> (R before F)", {
    #  R ---------->                        ------------> F  Primers
    # ====================================================== DNA Target
    #    ---------->************************------------>    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 480,
      forward_stop = 500,
      reverse_start = 100,
      reverse_stop = 120
    )

    expect_true(is.na(product_length))
  })

  it("returns NA with both primers like this: <--- (F before R)", {
    #  F <----------                        <------------ R  Primers
    # ====================================================== DNA Target
    #    <----------************************<------------    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 120,
      forward_stop = 100,
      reverse_start = 500,
      reverse_stop = 480
    )

    expect_true(is.na(product_length))
  })

  it("returns NA with both primers like this: <--- (F before R)", {
    #  R <----------                        <------------ F  Primers
    # ====================================================== DNA Target
    #    <----------************************<------------    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 500,
      forward_stop = 480,
      reverse_start = 120,
      reverse_stop = 100
    )

    expect_true(is.na(product_length))
  })

  it("returns NA with F <---, and R ---> (F before R)", {
    #  F <----------                        ------------> R  Primers
    # ====================================================== DNA Target
    #    <----------************************------------>    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 120,
      forward_stop = 100,
      reverse_start = 480,
      reverse_stop = 500
    )

    expect_true(is.na(product_length))
  })

  it("returns NA with F <---, and R ---> (R before F)", {
    #  R <----------                        ------------> F  Primers
    # ====================================================== DNA Target
    #    <----------************************------------>    Product
    #    ^         ^                        ^           ^
    #    100       120                      480         500
    product_length <- .calculate_product_length(
      forward_start = 480,
      forward_stop = 500,
      reverse_start = 120,
      reverse_stop = 100
    )

    expect_true(is.na(product_length))
  })
})

describe(".find_all_product_lengths", {
  it("joins forward and reverse hits and adds product_length column", {
    forward_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC002"),
      subject_gi = c("GI001", "GI002"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX002"),
      forward_start = c(100, 200),
      forward_stop = c(120, 220),
      forward_mismatch = c(0, 1)
    )

    reverse_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC002"),
      subject_gi = c("GI001", "GI002"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX002"),
      reverse_start = c(500, 700),
      reverse_stop = c(480, 680),
      reverse_mismatch = c(0, 2)
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_equal(nrow(result), 2)
    expect_equal(result$product_length, c(400, 500))
    expect_equal(result$forward_mismatch, c(0, 1))
    expect_equal(result$reverse_mismatch, c(0, 2))
  })

  it("handles many-to-many relationships (one forward to multiple reverse)", {
    forward_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      forward_start = 100,
      forward_stop = 120,
      forward_mismatch = 0
    )

    reverse_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC001", "ACC001"),
      subject_gi = c("GI001", "GI001", "GI001"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX001", "TAX001"),
      reverse_start = c(500, 600, 700),
      reverse_stop = c(480, 580, 680),
      reverse_mismatch = c(0, 1, 2)
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_equal(nrow(result), 3)
    expect_equal(result$product_length, c(400, 500, 600))
  })

  it("handles many-to-many relationships (multiple forward to one reverse)", {
    forward_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC001", "ACC001"),
      subject_gi = c("GI001", "GI001", "GI001"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX001", "TAX001"),
      forward_start = c(100, 200, 300),
      forward_stop = c(120, 220, 320),
      forward_mismatch = c(0, 1, 2)
    )

    reverse_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      reverse_start = 500,
      reverse_stop = 480,
      reverse_mismatch = 0
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_equal(nrow(result), 3)
    expect_equal(result$product_length, c(400, 300, 200))
  })

  it("handles many-to-many relationships (multiple forward to multiple reverse)", {
    forward_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC001"),
      subject_gi = c("GI001", "GI001"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX001"),
      forward_start = c(100, 200),
      forward_stop = c(120, 220),
      forward_mismatch = c(0, 1)
    )

    reverse_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC001"),
      subject_gi = c("GI001", "GI001"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX001"),
      reverse_start = c(500, 600),
      reverse_stop = c(480, 580),
      reverse_mismatch = c(0, 1)
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    # 2 forward x 2 reverse = 4 combinations
    expect_equal(nrow(result), 4)
    expect_equal(result$product_length, c(400, 500, 300, 400))
  })

  it("only includes matching records (inner join behavior)", {
    forward_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC002"),
      subject_gi = c("GI001", "GI002"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX002"),
      forward_start = c(100, 200),
      forward_stop = c(120, 220),
      forward_mismatch = c(0, 1)
    )

    reverse_hits <- data.frame(
      subject_accession_version = c("ACC001", "ACC003"),
      subject_gi = c("GI001", "GI003"),
      unique_subject_taxonomy_ids = c("TAX001", "TAX003"),
      reverse_start = c(500, 600),
      reverse_stop = c(480, 580),
      reverse_mismatch = c(0, 2)
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    # Only ACC001 matches
    expect_equal(nrow(result), 1)
    expect_equal(result$subject_accession_version, "ACC001")
    expect_equal(result$product_length, 400)
  })

  it("returns empty dataframe when no matches", {
    forward_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      forward_start = 100,
      forward_stop = 120,
      forward_mismatch = 0
    )

    reverse_hits <- data.frame(
      subject_accession_version = "ACC002",
      subject_gi = "GI002",
      unique_subject_taxonomy_ids = "TAX002",
      reverse_start = 500,
      reverse_stop = 480,
      reverse_mismatch = 0
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_equal(nrow(result), 0)
    expect_true("product_length" %in% colnames(result))
  })

  it("has the correct output columns", {
    forward_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      forward_start = 100,
      forward_stop = 120,
      forward_mismatch = 0
    )

    reverse_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      reverse_start = 500,
      reverse_stop = 480,
      reverse_mismatch = 2
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_true(all(
      c(
        "subject_accession_version",
        "subject_gi",
        "unique_subject_taxonomy_ids",
        "forward_start",
        "forward_stop",
        "forward_mismatch",
        "reverse_start",
        "reverse_stop",
        "reverse_mismatch",
        "product_length"
      ) %in%
        colnames(result)
    ))
  })

  it("includes NA product_length for invalid primer orientations", {
    forward_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      forward_start = 100,
      forward_stop = 120,
      forward_mismatch = 0
    )

    # Both primers pointing in same direction (invalid)
    reverse_hits <- data.frame(
      subject_accession_version = "ACC001",
      subject_gi = "GI001",
      unique_subject_taxonomy_ids = "TAX001",
      reverse_start = 480,
      reverse_stop = 500,
      reverse_mismatch = 0
    )

    result <- .find_all_product_lengths(forward_hits, reverse_hits)

    expect_equal(nrow(result), 1)
    expect_true(is.na(result$product_length))
  })
})

# TODO: add a couple integration tests for find_plausible_amplicon_coordinates
# that brings it all together
