describe("taxonomy_ids_to_taxonomy_table", {
  taxonomy_db <- testthat::test_path(
    "data",
    "small_test_db",
    "taxonomy.db"
  )

  it("returns rows of NA values for taxonomy IDs that don't exist in the DB", {
    taxonomy_ids <- c("1861737", "i_dont_exist", "1284384")

    result <- taxonomy_ids_to_taxonomy_table(
      taxonomy_ids,
      taxonomy_db,
      desired_taxa = c("domain", "family")
    )

    testthat::expect_equal(
      result,
      expected = tibble::tibble(
        taxonomy_id = taxonomy_ids,
        domain = c("Eukaryota", NA, "Eukaryota"),
        family = c("Stichococcaceae", NA, "Prasiolaceae")
      )
    )
  })
})
