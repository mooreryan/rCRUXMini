basic_blastdbcmd_result <- function() {
  # This data comes from the spreadsheet. Easier to deal with and update than
  # tibble/tribble directly.
  readr::read_tsv(
    I(
      "1	100	S101-1	w	1011	db_1
      2	100	S101-2	w	1012	db_1
      3	100	S101-3	w	1013	db_1
      4	200	S102-1	x	1021	db_1
      5	200	S102-2	x	1022	db_1
      6	300	S103-1	y	1031	db_1
      7	300	S103-2	y	1032	db_1
      8	400	S104-1	z	1041	db_2
      9	400	S104-2	z	1042	db_2
      10	400	S104-3	z	1043	db_2
      11	400	S104-4	z	1044	db_2"
    ),
    trim_ws = TRUE,
    col_names = c(
      "index",
      "blast_ordinal_id",
      "subject_accession_version",
      "sequence_hash_value",
      "unique_subject_taxonomy_ids",
      "blast_db_path"
    ),
    col_types = "iicccc"
  )
}

basic_blast_result <- function() {
  readr::read_tsv(
    I(
      "Q1	S101-1	1	1	0.0001	1	1	1	A	1011;1012;1013
      Q2	S102-1	2	2	0.0002	2	2	2	AA	1021;1022
      Q3	S102-1	3	3	0.0003	2	2	2	CC	1021;1022
      Q4	S103-1	4	4	0.0004	3	3	3	AAA	1031;1032
      Q5	S103-1	5	5	0.0005	3	3	3	CCC	1031;1032
      Q1	S103-1	6	6	0.0006	3	3	3	GGG	1031;1032
      Q2	S104-1	7	7	0.0007	4	4	4	AAAA	1041;1042;1043;1044
      Q3	S104-1	8	8	0.0008	4	4	4	CCCC	1041;1042;1043;1044
      Q4	S104-1	9	9	0.0009	4	4	4	GGGG	1041;1042;1043;1044
      Q5	S104-1	10	10	0.001	4	4	4	TTTT	1041;1042;1043;1044"
    ),
    trim_ws = TRUE,
    col_names = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids"
    ),
    col_types = c("ccdidiiicc")
  )
}

# The rust parsers spit out blast results with degapped columns
basic_with_degapped_blast_result <- function() {
  readr::read_tsv(
    I(
      "Q1	S101-1	1	1	0.0001	1	1	1	A	1011;1012;1013	1	A
      Q2	S102-1	2	2	0.0002	2	2	2	AA	1021;1022	2	AA
      Q3	S102-1	3	3	0.0003	2	2	2	CC	1021;1022	2	CC
      Q4	S103-1	4	4	0.0004	3	3	3	AAA	1031;1032	3	AAA
      Q5	S103-1	5	5	0.0005	3	3	3	CCC	1031;1032	3	CCC
      Q1	S103-1	6	6	0.0006	3	3	3	GGG	1031;1032	3	GGG
      Q2	S104-1	7	7	0.0007	4	4	4	AAAA	1041;1042;1043;1044	4	AAAA
      Q3	S104-1	8	8	0.0008	4	4	4	CCCC	1041;1042;1043;1044	4	CCCC
      Q4	S104-1	9	9	0.0009	4	4	4	GGGG	1041;1042;1043;1044	4	GGGG
      Q5	S104-1	10	10	0.001	4	4	4	TTTT	1041;1042;1043;1044	4	TTTT"
    ),
    trim_ws = TRUE,
    col_names = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids",
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    ),
    col_types = c("ccdidiiicccc")
  )
}


expanded_result <- function() {
  readr::read_tsv(
    I(
      "Q1	S101-1	1	1	0.0001	1	1	1	A	1011
      Q1	S101-2	1	1	0.0001	1	1	1	A	1012
      Q1	S101-3	1	1	0.0001	1	1	1	A	1013
      Q2	S102-1	2	2	0.0002	2	2	2	AA	1021
      Q2	S102-2	2	2	0.0002	2	2	2	AA	1022
      Q3	S102-1	3	3	0.0003	2	2	2	CC	1021
      Q3	S102-2	3	3	0.0003	2	2	2	CC	1022
      Q4	S103-1	4	4	0.0004	3	3	3	AAA	1031
      Q4	S103-2	4	4	0.0004	3	3	3	AAA	1032
      Q5	S103-1	5	5	0.0005	3	3	3	CCC	1031
      Q5	S103-2	5	5	0.0005	3	3	3	CCC	1032
      Q1	S103-1	6	6	0.0006	3	3	3	GGG	1031
      Q1	S103-2	6	6	0.0006	3	3	3	GGG	1032
      Q2	S104-1	7	7	0.0007	4	4	4	AAAA	1041
      Q2	S104-2	7	7	0.0007	4	4	4	AAAA	1042
      Q2	S104-3	7	7	0.0007	4	4	4	AAAA	1043
      Q2	S104-4	7	7	0.0007	4	4	4	AAAA	1044
      Q3	S104-1	8	8	0.0008	4	4	4	CCCC	1041
      Q3	S104-2	8	8	0.0008	4	4	4	CCCC	1042
      Q3	S104-3	8	8	0.0008	4	4	4	CCCC	1043
      Q3	S104-4	8	8	0.0008	4	4	4	CCCC	1044
      Q4	S104-1	9	9	0.0009	4	4	4	GGGG	1041
      Q4	S104-2	9	9	0.0009	4	4	4	GGGG	1042
      Q4	S104-3	9	9	0.0009	4	4	4	GGGG	1043
      Q4	S104-4	9	9	0.0009	4	4	4	GGGG	1044
      Q5	S104-1	10	10	0.001	4	4	4	TTTT	1041
      Q5	S104-2	10	10	0.001	4	4	4	TTTT	1042
      Q5	S104-3	10	10	0.001	4	4	4	TTTT	1043
      Q5	S104-4	10	10	0.001	4	4	4	TTTT	1044"
    ),
    trim_ws = TRUE,
    col_names = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids"
    ),
    col_types = c("ccdidiiicc")
  )
}

expanded_with_degapped_result <- function() {
  readr::read_tsv(
    I(
      "Q1	S101-1	1	1	0.0001	1	1	1	A	1011	1	A
      Q1	S101-2	1	1	0.0001	1	1	1	A	1012	1	A
      Q1	S101-3	1	1	0.0001	1	1	1	A	1013	1	A
      Q2	S102-1	2	2	0.0002	2	2	2	AA	1021	2	AA
      Q2	S102-2	2	2	0.0002	2	2	2	AA	1022	2	AA
      Q3	S102-1	3	3	0.0003	2	2	2	CC	1021	2	CC
      Q3	S102-2	3	3	0.0003	2	2	2	CC	1022	2	CC
      Q4	S103-1	4	4	0.0004	3	3	3	AAA	1031	3	AAA
      Q4	S103-2	4	4	0.0004	3	3	3	AAA	1032	3	AAA
      Q5	S103-1	5	5	0.0005	3	3	3	CCC	1031	3	CCC
      Q5	S103-2	5	5	0.0005	3	3	3	CCC	1032	3	CCC
      Q1	S103-1	6	6	0.0006	3	3	3	GGG	1031	3	GGG
      Q1	S103-2	6	6	0.0006	3	3	3	GGG	1032	3	GGG
      Q2	S104-1	7	7	0.0007	4	4	4	AAAA	1041	4	AAAA
      Q2	S104-2	7	7	0.0007	4	4	4	AAAA	1042	4	AAAA
      Q2	S104-3	7	7	0.0007	4	4	4	AAAA	1043	4	AAAA
      Q2	S104-4	7	7	0.0007	4	4	4	AAAA	1044	4	AAAA
      Q3	S104-1	8	8	0.0008	4	4	4	CCCC	1041	4	CCCC
      Q3	S104-2	8	8	0.0008	4	4	4	CCCC	1042	4	CCCC
      Q3	S104-3	8	8	0.0008	4	4	4	CCCC	1043	4	CCCC
      Q3	S104-4	8	8	0.0008	4	4	4	CCCC	1044	4	CCCC
      Q4	S104-1	9	9	0.0009	4	4	4	GGGG	1041	4	GGGG
      Q4	S104-2	9	9	0.0009	4	4	4	GGGG	1042	4	GGGG
      Q4	S104-3	9	9	0.0009	4	4	4	GGGG	1043	4	GGGG
      Q4	S104-4	9	9	0.0009	4	4	4	GGGG	1044	4	GGGG
      Q5	S104-1	10	10	0.001	4	4	4	TTTT	1041	4	TTTT
      Q5	S104-2	10	10	0.001	4	4	4	TTTT	1042	4	TTTT
      Q5	S104-3	10	10	0.001	4	4	4	TTTT	1043	4	TTTT
      Q5	S104-4	10	10	0.001	4	4	4	TTTT	1044	4	TTTT"
    ),
    trim_ws = TRUE,
    col_names = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids",
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    ),
    col_types = c("ccdidiiicccc")
  )
}

test_that(".add_reference_subject_accession_versions works", {
  blastdbcmd_result <- basic_blastdbcmd_result()

  expected <- c(
    rep("S101-1", 3),
    rep("S102-1", 2),
    rep("S103-1", 2),
    rep("S104-1", 4)
  )

  actual <- blastdbcmd_result |>
    .add_reference_subject_accession_versions() |>
    dplyr::pull(.data$reference_subject_accession_version)

  expect_equal(actual, expected)
})

# devtools::test_active_file("tests/testthat/test-expand_multi_tax_ids.R")

test_that(".expand_hits works", {
  blastdbcmd_result <- basic_blastdbcmd_result() |>
    .add_reference_subject_accession_versions()
  blast_result <- basic_blast_result()

  result <- .expand_hits(
    blastdbcmd_result = blastdbcmd_result,
    blast_result = blast_result
  ) |>
    dplyr::arrange(.data$subject_accession_version, .data$query_accession)

  expected <- expanded_result() |>
    dplyr::arrange(.data$subject_accession_version, .data$query_accession)

  expect_equal(object = result, expected = expected)
})


describe("expand_multi_tax_ids2", {
  blast_db_paths <- c("db_1", "db_2")

  it("works with sequential plan", {
    # Ensure the sequential plan is set for this test
    with(future::plan(future::sequential), local = TRUE)

    result <- expand_multi_tax_ids2(
      blastdbcmd = testthat::test_path("scripts/blastdbcmd"),
      hits_with_multiple_taxonomy_ids = basic_with_degapped_blast_result(),
      blast_db_paths = blast_db_paths
    ) |>
      dplyr::arrange(.data$subject_accession_version, .data$query_accession)

    expected <- expanded_with_degapped_result() |>
      dplyr::arrange(.data$subject_accession_version, .data$query_accession)

    expect_equal(result, expected = expected)
  })

  it("works with multisession plan", {
    # Ensure the multisession plan is set for this test
    with(future::plan(future::multisession, workers = 3), local = TRUE)

    result <- expand_multi_tax_ids2(
      blastdbcmd = testthat::test_path("scripts/blastdbcmd"),
      hits_with_multiple_taxonomy_ids = basic_with_degapped_blast_result(),
      blast_db_paths = blast_db_paths
    ) |>
      dplyr::arrange(.data$subject_accession_version, .data$query_accession)

    expected <- expanded_with_degapped_result() |>
      dplyr::arrange(.data$subject_accession_version, .data$query_accession)

    expect_equal(result, expected = expected)
  })

  it("handles flaky blastdbcmd", {
    # Ensure the multisession plan is set for this test
    with(future::plan(future::multisession, workers = 3), local = TRUE)

    testthat::expect_error(
      {
        expand_multi_tax_ids2(
          blastdbcmd = testthat::test_path("scripts/blastdbcmd"),
          hits_with_multiple_taxonomy_ids = basic_with_degapped_blast_result(),
          blast_db_paths = c(blast_db_paths, "db_3")
        )
      },
      class = "rcrux_mini_error"
    )
  })
})
