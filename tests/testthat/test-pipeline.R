library(testthat)
library(checkmate)
library(dplyr)
library(readr)

# Pass in the unique ID and it returns a function that that will call
# expect_snapshot_file, but with a unique ID appended to the name of the file.
# Use this to keep the file snapshots unique to the individual test even with a
# single test file.
expect_snapshot_file_factory <- function(id) {
  function(path) {
    name <- stringr::str_glue("{id}__{basename(path)}")
    testthat::expect_snapshot_file(path = path, name = name)
  }
}

test_that("to_fasta_string given a single sequence", {
  # Arrange
  input_sequence <- "AAAA"
  expected <- ">xyz_1\nAAAA"

  # Act
  result <- to_fasta_string("xyz", input_sequence)

  # Assert
  testthat::expect_equal(result, expected)
})

test_that("to_fasta_string given a multiple sequences", {
  # Arrange
  input_sequence <- c("AAAA", "TTTT")
  expected <- ">xyz_1\nAAAA\n>xyz_2\nTTTT"

  # Act
  result <- to_fasta_string("xyz", input_sequence)

  # Assert
  testthat::expect_equal(result, expected)
})

test_that("parse_primer_blast_results given basic valid input returns parsed data frame", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi2"),
    saccver = c("ACC001", "ACC001"),
    mismatch = c(1, 2),
    sstart = c(100, 200),
    send = c(120, 220),
    staxids = c("12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  testthat::expect_equal(nrow(result), 2)
  testthat::expect_true(all(
    c("primer_forward", "primer_reverse") %in% result$qseqid
  ))

  unlink(input_file)
})

test_that("parse_primer_blast_results given empty file throws error", {
  # Arrange
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(data.frame(), input_file)

  # Act & Assert
  testthat::expect_error(
    parse_primer_blast_results(input_file),
    "NO PRIMER BLAST RESULTS"
  )

  unlink(input_file)
})

test_that("parse_primer_blast_results given sequences without both primers filters out incomplete pairs", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi2", "gi3"),
    saccver = c("ACC001", "ACC002", "ACC001"), # ACC002 only has forward
    mismatch = c(1, 1, 2),
    sstart = c(100, 100, 200),
    send = c(120, 120, 220),
    staxids = c("12345", "12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_accession <- "ACC001"

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  testthat::expect_true(all(result$saccver == expected_accession))
  testthat::expect_false("ACC002" %in% result$saccver)

  unlink(input_file)
})

test_that("parse_primer_blast_results given high mismatch sequences filters out poor matches", {
  # Arrange
  input_data <- data.frame(
    qseqid = c(
      "primer_forward",
      "primer_reverse",
      "primer_forward",
      "primer_reverse"
    ),
    sgi = c("gi1", "gi2", "gi3", "gi4"),
    saccver = c("ACC001", "ACC001", "ACC002", "ACC002"),
    mismatch = c(1, 2, 5, 4), # ACC002 has mismatches >= 4
    sstart = c(100, 200, 100, 200),
    send = c(120, 220, 120, 220),
    staxids = c("12345", "12345", "12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_accession <- "ACC001"
  max_allowed_mismatches <- 4

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  testthat::expect_true(all(result$saccver == expected_accession))
  testthat::expect_true(all(result$mismatch < max_allowed_mismatches))

  unlink(input_file)
})

test_that("parse_primer_blast_results given multiple matches per primer selects best by minimum mismatches", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi2", "gi3"),
    saccver = c("ACC001", "ACC001", "ACC001"),
    mismatch = c(3, 1, 2), # Two forwards with different mismatches
    sstart = c(100, 100, 200), # Same start position for forwards
    send = c(120, 120, 220),
    staxids = c("12345", "12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_best_mismatch <- 1

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  forward_result <- result[result$qseqid == "primer_forward", ]
  testthat::expect_equal(forward_result$mismatch, expected_best_mismatch)

  unlink(input_file)
})

test_that("parse_primer_blast_results given duplicate start positions removes duplicates", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi2", "gi3"),
    saccver = c("ACC001", "ACC001", "ACC001"),
    mismatch = c(1, 1, 2), # Same mismatches
    sstart = c(100, 100, 200), # Duplicate start positions
    send = c(120, 125, 220), # Different end positions
    staxids = c("12345", "12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_forward_count <- 1

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  forward_results <- result[result$qseqid == "primer_forward", ]
  testthat::expect_equal(nrow(forward_results), expected_forward_count)

  unlink(input_file)
})

test_that("parse_primer_blast_results given duplicate end positions removes duplicates", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse", "primer_reverse"),
    sgi = c("gi1", "gi2", "gi3"),
    saccver = c("ACC001", "ACC001", "ACC001"),
    mismatch = c(1, 2, 2), # Same mismatches for reverses
    sstart = c(100, 200, 205), # Different start positions for reverses
    send = c(120, 220, 220), # Same end position for reverses
    staxids = c("12345", "12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_reverse_count <- 1

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  reverse_results <- result[result$qseqid == "primer_reverse", ]
  testthat::expect_equal(nrow(reverse_results), expected_reverse_count)

  unlink(input_file)
})

test_that("parse_primer_blast_results given multiple accessions returns results for all valid accessions", {
  # Arrange
  input_data <- data.frame(
    qseqid = c(
      "primer_forward",
      "primer_reverse",
      "primer_forward",
      "primer_reverse"
    ),
    sgi = c("gi1", "gi2", "gi3", "gi4"),
    saccver = c("ACC001", "ACC001", "ACC002", "ACC002"),
    mismatch = c(1, 2, 0, 3),
    sstart = c(100, 200, 150, 250),
    send = c(120, 220, 170, 270),
    staxids = c("12345", "12345", "67890", "67890")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_accessions <- c("ACC001", "ACC002")
  expected_row_count <- 4

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  testthat::expect_true(all(expected_accessions %in% result$saccver))
  testthat::expect_equal(nrow(result), expected_row_count)

  unlink(input_file)
})

test_that("parse_primer_blast_results given valid input maintains expected column structure", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi2"),
    saccver = c("ACC001", "ACC001"),
    mismatch = c(1, 2),
    sstart = c(100, 200),
    send = c(120, 220),
    staxids = c("12345", "12345")
  )
  input_file <- tempfile(fileext = ".tsv")
  write_tsv(input_data, input_file, col_names = FALSE)
  expected_columns <- c(
    "qseqid",
    "sgi",
    "saccver",
    "mismatch",
    "sstart",
    "send",
    "staxids"
  )

  # Act
  result <- parse_primer_blast_results(input_file)

  # Assert
  testthat::expect_true(all(expected_columns %in% names(result)))

  unlink(input_file)
})


##############################################################################
# find_plausible_amplicon_coordinates ###################################################
##############################################################################

test_that("find_plausible_amplicon_coordinates validates input parameters correctly", {
  # Create minimal valid data frame
  valid_data <- data.frame(
    qseqid = c("forward_1", "reverse_1"),
    sgi = c("gi1", "gi2"),
    saccver = c("ACC001", "ACC001"),
    mismatch = c(1, 2),
    sstart = c(100, 200),
    send = c(120, 220),
    staxids = c("12345", "12345")
  )

  # Test invalid data frame structure
  expect_error(
    find_plausible_amplicon_coordinates(
      data.frame(),
      minimum_length = 50,
      maximum_length = 500
    ),
    "have at least 1 row"
  )

  # Test missing required columns
  expect_error(
    find_plausible_amplicon_coordinates(
      valid_data[, -1], # remove qseqid column
      minimum_length = 50,
      maximum_length = 500
    ),
    "[Nn]ames must include"
  )

  # # Test invalid length parameters
  expect_error(
    find_plausible_amplicon_coordinates(
      valid_data,
      minimum_length = 0, # should be >= 1
      maximum_length = 500
    ),
    "minimum_length"
  )

  expect_error(
    find_plausible_amplicon_coordinates(
      valid_data,
      # min is > max
      minimum_length = 500,
      maximum_length = 100
    ),
    "minimum_length"
  )
})

test_that("find_plausible_amplicon_coordinates calculates product length correctly for standard orientation", {
  # Arrange: Forward primer before reverse primer
  # F ---------->                        <------------ R   Primers
  # ====================================================== DNA Target
  # Expected amplicon length = reverse_start - forward_start = 200 - 100 = 100
  input_data <- data.frame(
    qseqid = c("forward_1", "reverse_1"),
    sgi = c("gi1", "gi1"),
    saccver = c("ACC001", "ACC001"),
    mismatch = c(1, 2),
    sstart = c(100, 200), # forward starts at 100, reverse starts at 200
    send = c(120, 180), # forward ends at 120, reverse ends at 180
    staxids = c("12345", "12345")
  )
  expected_length <- 100

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$product_length, expected_length)
  expect_equal(result$accession, "ACC001")
})

test_that("find_plausible_amplicon_coordinates calculates product length correctly for reverse orientation", {
  # Arrange: Reverse primer before forward primer
  # R ---------->                        <------------ F   Primers
  # ====================================================== DNA Target
  # Expected amplicon length = forward_start - reverse_start = 200 - 100 = 100
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = c("gi1", "gi1"),
    saccver = c("ACC001", "ACC001"),
    mismatch = c(1, 2),
    sstart = c(200, 100), # forward starts at 200, reverse starts at 100
    send = c(180, 120), # forward ends at 180, reverse ends at 120
    staxids = c("12345", "12345")
  )
  expected_length <- 100

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$product_length, expected_length)
})

test_that("find_plausible_amplicon_coordinates filters amplicons by length correctly", {
  # Arrange: Create amplicons of different lengths
  input_data <- data.frame(
    qseqid = c(
      "primer_forward",
      "primer_reverse",
      "primer_forward",
      "primer_reverse"
    ),
    sgi = c("gi1", "gi1", "gi2", "gi2"),
    saccver = c("ACC001", "ACC001", "ACC002", "ACC002"),
    mismatch = c(1, 2, 1, 2),
    sstart = c(100, 200, 100, 400), # ACC001: 200-100=100, ACC002: 400-100=300
    send = c(120, 180, 120, 380),
    staxids = c("12345", "12345", "67890", "67890")
  )

  # Act: Set length filter to exclude the 300bp amplicon
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 299 # Should exclude ACC002
  )

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$accession, "ACC001")
  expect_equal(result$product_length, 100)

  # Act: Set length filter to exclude the 100bp amplicon
  result <- find_plausible_amplicon_coordinates(
    input_data,
    # Should exclude ACC002
    minimum_length = 101,
    maximum_length = 301
  )

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$accession, "ACC002")
  expect_equal(result$product_length, 300)
})

test_that("find_plausible_amplicon_coordinates selects best primer combinations per accession", {
  # Arrange: Multiple forward and reverse hits for same accession
  input_data <- data.frame(
    qseqid = c(
      "forward_1",
      "reverse_1",
      # Names don't matter
      "forward_2",
      "reverse_2",
      # This one would make the biggest amplicon, but it is in the wrong orientation
      "reverse_3"
    ),
    sgi = "gi1",
    saccver = "ACC001",
    mismatch = c(1, 2, 1, 2, 1),
    sstart = c(100, 150, 200, 250, 500),
    send = c(120, 170, 180, 230, 525),
    staxids = "12345"
  )

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert: Should keep the longest combination in the correct orientation
  expect_equal(nrow(result), 1)
  expect_equal(result$accession, "ACC001")
  expect_true(result$product_length == 150)
})

test_that("find_plausible_amplicon_coordinates returns empty result when no valid amplicons found", {
  # Arrange: Primers that don't meet length criteria
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = "gi1",
    saccver = "ACC001",
    mismatch = c(1, 2),
    sstart = c(100, 110), # Very short amplicon length = 10
    send = c(105, 105),
    staxids = "12345"
  )

  # Act: Set minimum length higher than possible amplicon
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_equal(nrow(result), 0)
})

test_that("find_plausible_amplicon_coordinates overlapping primers in the correct orientation and size yield plausible amplicons", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = "gi1",
    saccver = "ACC001",
    mismatch = c(1, 2),
    sstart = c(100, 200),
    send = c(250, 180), # forward_stop > reverse_start
    staxids = "12345"
  )

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_equal(nrow(result), 1)
})

test_that("find_plausible_amplicon_coordinates handles multiple accessions correctly", {
  # Arrange: Multiple accessions with valid amplicons and some with no amplicons
  input_data <- data.frame(
    qseqid = c(
      # ACC001
      "primer_forward",
      "primer_reverse",
      # ACC002
      "primer_forward",
      "primer_reverse",
      # ACC003 has only forward primer
      "foward_1",
      # ACC004 has only reverse primer
      "reverse_1"
    ),
    sgi = c("gi1", "gi1", "gi2", "gi2", "gi3", "gi4"),
    saccver = c("ACC001", "ACC001", "ACC002", "ACC002", "ACC003", "ACC004"),
    mismatch = c(1, 2, 0, 1, 0, 1),
    sstart = c(100, 200, 150, 250, 100, 200),
    send = c(120, 180, 170, 230, 120, 180),
    staxids = c("12345", "12345", "67890", "67890", "3456", "4567")
  )

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_equal(nrow(result), 2)
  expect_true(all(result$accession %in% c("ACC001", "ACC002")))
  expect_equal(result$product_length, c(100, 100))
})

test_that("find_plausible_amplicon_coordinates has required columns in output", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = "gi1",
    saccver = "ACC001",
    mismatch = c(1, 2),
    sstart = c(100, 200),
    send = c(120, 180),
    staxids = "12345"
  )
  expected_columns <- c(
    "accession",
    "gi",
    "staxids",
    "forward_start",
    "forward_stop",
    "forward_mismatch",
    "reverse_start",
    "reverse_stop",
    "reverse_mismatch",
    "product_length"
  )

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert
  expect_names(names(result), must.include = expected_columns)
  expect_equal(result$accession, "ACC001")
  expect_equal(result$gi, "gi1")
  expect_equal(result$staxids, "12345")
})

test_that("find_plausible_amplicon_coordinates handles missing forward or reverse primers", {
  # Arrange: Only forward primers, no reverse
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_forward"),
    sgi = c("gi1", "gi2"),
    saccver = c("ACC001", "ACC002"),
    mismatch = c(1, 2),
    sstart = c(100, 150),
    send = c(120, 170),
    staxids = c("12345", "67890")
  )

  # Act
  result <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500
  )

  # Assert: Should return empty because amplicon requires both forward and
  # reverse primers
  expect_equal(nrow(result), 0)
})

test_that("find_plausible_amplicon_coordinates respects maximum_mismatches parameter", {
  # Arrange
  input_data <- data.frame(
    qseqid = c("primer_forward", "primer_reverse"),
    sgi = "gi1",
    saccver = c("ACC001", "ACC001"),
    mismatch = c(3, 3), # exactly at boundary
    sstart = c(100, 200),
    send = c(120, 180),
    staxids = "12345"
  )

  # Act with different maximum_mismatches values
  result_strict <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500,
    maximum_mismatches = 2 # Should exclude data
  )

  result_permissive <- find_plausible_amplicon_coordinates(
    input_data,
    minimum_length = 50,
    maximum_length = 500,
    maximum_mismatches = 4 # Should include data
  )

  # Assert
  expect_equal(nrow(result_strict), 0)
  expect_equal(nrow(result_permissive), 1)
})

##############################################################################
# accession_to_taxonomy ######################################################
##############################################################################

test_that("accession_to_taxonomy validates input parameters correctly", {
  # Test empty data frame
  expect_error(
    accession_to_taxonomy(
      data.frame(),
      taxonomy_db_path = file.path(
        test_path(),
        "data",
        "small_test_db",
        "taxonomy.db"
      )
    ),
    "have at least 1 row"
  )

  # Test missing accession column
  expect_error(
    accession_to_taxonomy(
      data.frame(other_col = "value"),
      taxonomy_db_path = file.path(
        test_path(),
        "data",
        "small_test_db",
        "taxonomy.db"
      )
    ),
    "[Nn]ames must include.*accession"
  )

  # Test non-existent taxonomy database
  expect_error(
    accession_to_taxonomy(
      data.frame(accession = "test"),
      taxonomy_db_path = "/non/existent/path.db"
    ),
    "[Ff]ile does not exist"
  )
})

test_that("accession_to_taxonomy handles valid accessions correctly", {
  # Arrange
  input_data <- data.frame(
    accession = c("sequence_0", "sequence_1"),
    other_col = c("value1", "value2")
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_equal(nrow(result), 2)
  expect_true(all(
    c("accession", "other_col", "taxonomy_id") %in% names(result)
  ))

  # Check that valid accessions have non-NA taxonomy data
  sequence_0_row <- result[result$accession == "sequence_0", ]
  sequence_1_row <- result[result$accession == "sequence_1", ]

  expect_false(is.na(sequence_0_row$phylum))
  expect_false(is.na(sequence_1_row$phylum))
  expect_false(is.na(sequence_0_row$taxonomy_id))
  expect_false(is.na(sequence_1_row$taxonomy_id))
})

test_that("accession_to_taxonomy handles invalid accessions correctly", {
  # Arrange
  input_data <- data.frame(
    accession = c("invalid_accession_123", "another_invalid_456"),
    other_col = c("value1", "value2")
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_equal(nrow(result), 2)

  # Invalid accessions should have NA taxonomy data
  expect_true(all(is.na(result$phylum)))
  expect_true(all(is.na(result$taxonomy_id)))
})

test_that("accession_to_taxonomy handles mixed valid and invalid accessions", {
  # Arrange
  input_data <- data.frame(
    accession = c("sequence_0", "invalid_accession", "sequence_1"),
    other_col = c("value1", "value2", "value3")
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_equal(nrow(result), 3)

  # Valid accessions should have taxonomy data
  valid_rows <- result[result$accession %in% c("sequence_0", "sequence_1"), ]
  expect_true(all(!is.na(valid_rows$phylum)))
  expect_true(all(!is.na(valid_rows$taxonomy_id)))

  # Invalid accession should have NA taxonomy data
  invalid_row <- result[result$accession == "invalid_accession", ]
  expect_true(is.na(invalid_row$phylum))
  expect_true(is.na(invalid_row$taxonomy_id))
})

test_that("accession_to_taxonomy preserves original columns and adds taxonomy columns", {
  # Arrange
  input_data <- data.frame(
    accession = "sequence_0",
    gi = "gi123",
    other_data = "test_value",
    numeric_col = 1
  )

  expected_original_cols <- c("accession", "gi", "other_data", "numeric_col")
  expected_taxonomy_cols <- c(
    "taxonomy_id",
    "species",
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "family",
    "subfamily",
    "genus",
    "infraorder",
    "subcohort",
    "superorder",
    "superfamily",
    "tribe",
    "subspecies",
    "subgenus",
    "species.group",
    "parvorder",
    "varietas"
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_true(all(expected_original_cols %in% names(result)))
  expect_true(all(expected_taxonomy_cols %in% names(result)))

  # Original data should be preserved
  expect_equal(result$accession, "sequence_0")
  expect_equal(result$gi, "gi123")
  expect_equal(result$other_data, "test_value")
  expect_equal(result$numeric_col, 1)
})

test_that("accession_to_taxonomy preserves row count with original input", {
  # Arrange
  input_data <- data.frame(
    accession = c("sequence_0", "invalid_123", "sequence_1", "another_invalid"),
    gi = c("gi1", "gi2", "gi3", "gi4"),
    product_length = c(100, 200, 300, 400)
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_equal(nrow(result), nrow(input_data))
  expect_equal(sort(result$accession), sort(input_data$accession))
})

test_that("accession_to_taxonomy works with single accession", {
  # Arrange
  input_data <- tibble::tibble(
    accession = "sequence_0",
    gi = "gi123"
  )

  # Act
  result <- accession_to_taxonomy(
    input_data,
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    )
  )

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$accession, "sequence_0")
  expect_false(is.na(result$taxonomy_id))
})

##############################################################################
# distinct_taxonomic_ranks ###################################################
##############################################################################

test_that("distinct_taxonomic_ranks validates input parameters correctly", {
  # Test empty data frame
  expect_error(
    distinct_taxonomic_ranks(data.frame()),
    "have at least 1 row"
  )

  # Test missing required taxonomic columns
  incomplete_data <- data.frame(
    accession = c("ACC001", "ACC002"),
    superkingdom = c("Bacteria", "Bacteria"),
    phylum = c("Proteobacteria", "Firmicutes")
    # Missing class, order, family, genus, species
  )

  expect_error(
    distinct_taxonomic_ranks(incomplete_data),
    "[Nn]ames must include"
  )
})

test_that("distinct_taxonomic_ranks counts distinct values correctly for single row", {
  # Arrange
  input_data <- data.frame(
    accession = "ACC001",
    superkingdom = "SK1",
    phylum = "P1",
    class = "C1",
    order = "O1",
    family = "F1",
    genus = "G1",
    species = "S1"
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$superkingdom, 1)
  expect_equal(result$phylum, 1)
  expect_equal(result$class, 1)
  expect_equal(result$order, 1)
  expect_equal(result$family, 1)
  expect_equal(result$genus, 1)
  expect_equal(result$species, 1)
})

test_that("distinct_taxonomic_ranks counts distinct values correctly for multiple identical rows", {
  # Arrange - multiple rows with identical taxonomy
  input_data <- tibble::tibble(
    accession = c("ACC001", "ACC002", "ACC003"),
    superkingdom = "SK1",
    phylum = "P1",
    class = "C1",
    order = "O1",
    family = "F1",
    genus = "G1",
    species = "S1"
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert - should count as 1 distinct value for each level
  expect_equal(nrow(result), 1)
  expect_equal(result$superkingdom, 1)
  expect_equal(result$phylum, 1)
  expect_equal(result$class, 1)
  expect_equal(result$order, 1)
  expect_equal(result$family, 1)
  expect_equal(result$genus, 1)
  expect_equal(result$species, 1)
})

test_that("distinct_taxonomic_ranks counts distinct values correctly for diverse taxonomy", {
  # Arrange - multiple rows with different taxonomy at various levels
  input_data <- data.frame(
    accession = c("ACC001", "ACC002", "ACC003", "ACC004"),
    superkingdom = c("SK1", "SK1", "SK2", "SK1"),
    phylum = c("P1", "P2", "P3", "P4"),
    class = c("C1", "C2", "C3", "C1"),
    order = c("O1", "O2", "O3", "O1"),
    family = c("F1", "F2", "F3", "F1"),
    genus = c("G1", "G2", "G3", "G1"),
    species = c("S1", "S2", "S3", "S4")
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Because all phylums are different, but some of the lower orders are the
  # same, this test will show that this function doesn't match how the other
  # LCA script works.

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(result$superkingdom, 2)
  expect_equal(result$phylum, 4)
  expect_equal(result$class, 3)
  expect_equal(result$order, 3)
  expect_equal(result$family, 3)
  expect_equal(result$genus, 3)
  expect_equal(result$species, 4)
})

test_that("distinct_taxonomic_ranks handles NA values correctly", {
  # Arrange - data with some NA values in taxonomy
  input_data <- data.frame(
    accession = c("ACC001", "ACC002", "ACC003"),
    superkingdom = c("SK1", "SK1", NA),
    phylum = c("P1", NA, "P2"),
    class = c(NA, "C1", "C2"),
    order = c("O1", "O2", NA),
    family = c("F1", NA, NA),
    genus = c(NA, "G1", "G2"),
    species = c("S1", NA, "S2")
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert: NA values should be counted as distinct values
  expect_equal(nrow(result), 1)
  expect_equal(result$superkingdom, 2)
  expect_equal(result$phylum, 3)
  expect_equal(result$class, 3)
  expect_equal(result$order, 3)
  expect_equal(result$family, 2)
  expect_equal(result$genus, 3)
  expect_equal(result$species, 3)
})

test_that("distinct_taxonomic_ranks preserves only taxonomic level columns in output", {
  # Arrange - data with extra columns that should not appear in output
  input_data <- data.frame(
    accession = c("ACC001", "ACC002"),
    gi = c("gi1", "gi2"),
    product_length = c(100, 200),
    forward_start = c(50, 75),
    superkingdom = c("SK1", "SK1"),
    phylum = c("P1", "P2"),
    class = c("C1", "C2"),
    order = c("O1", "O2"),
    family = c("F1", "F2"),
    genus = c("G1", "G2"),
    species = c("S1", "S2")
  )

  expected_columns <- c(
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert
  expect_names(sort(names(result)), identical.to = sort(expected_columns))
})

test_that("distinct_taxonomic_ranks returns correct data types", {
  # Arrange
  input_data <- data.frame(
    superkingdom = c("S1", "S2"),
    phylum = c("P1", "P2"),
    class = c("C1", "C2"),
    order = c("O1", "O2"),
    family = c("F1", "F2"),
    genus = c("G1", "G2"),
    species = c("S1", "S2")
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert - all columns should be numeric (counts)
  expect_true(all(sapply(result, is.numeric)))
  expect_true(all(sapply(result, function(x) all(x >= 0))))
})

test_that("distinct_taxonomic_ranks works with minimal required data", {
  # Arrange - data frame with only the required taxonomic columns
  input_data <- data.frame(
    superkingdom = "SK1",
    phylum = "P1",
    class = "C1",
    order = "O1",
    family = "F1",
    genus = "G1",
    species = "S1"
  )

  # Act
  result <- distinct_taxonomic_ranks(input_data)

  # Assert
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 7)
  expect_true(all(result == 1))
})

##############################################################################
# get_seeds ##################################################################
##############################################################################

test_that("get_seeds", {
  expect_snapshot <- expect_snapshot_file_factory(
    "aa03fc3b-585b-4bb5-8cc9-0fda30bfc3ba"
  )

  output_directory_path_top <- tempfile()
  output_directory_path <- file.path(
    output_directory_path_top,
    "v2-get_seeds-out"
  )

  result <- get_seeds(
    forward_primers = c(
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC",
      "NN"
    ),
    reverse_primers = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT",
    output_directory_path = output_directory_path,
    blast_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "generated_sequences"
    ),
    taxonomy_db_path = file.path(
      test_path(),
      "data",
      "small_test_db",
      "taxonomy.db"
    ),
    random_seed = 254235
  )

  expect_snapshot(file.path(output_directory_path, "primers.fasta"))
  expect_snapshot(file.path(output_directory_path, "primer_blast.tsv"))

  print(system2(
    "ls",
    args = c("-l", output_directory_path),
    stdout = TRUE,
    stderr = TRUE
  ))

  print(result)
})
