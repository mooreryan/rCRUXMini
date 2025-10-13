# NOTE: the result is in the same order as the input taxonomy ids, and includes
# the taxonomy ID for each row
taxonomy_ids_to_taxonomy_table <- function(
  taxonomy_ids,
  taxonomy_db_path,
  desired_taxa = c(
    # This is the order from the original rCRUX
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
    "species group",
    "parvorder",
    "varietas"
  )
) {
  taxonomy_matrix <- taxonomizr_get_taxonomy(
    taxonomy_ids = taxonomy_ids,
    taxonomy_db_path = taxonomy_db_path,
    desired_taxa = desired_taxa
  )
  checkmate::assert_matrix(taxonomy_matrix, nrows = length(taxonomy_ids))

  result <- tibble::tibble(
    taxonomy_id = taxonomy_ids,
    data.frame(taxonomy_matrix)
  )
  checkmate::assert_data_frame(result, nrows = length(taxonomy_ids))

  result
}

# NOTE: If a taxonomy ID is present in the taxonomy_ids input, but it is _not_
# present in the taxonomy DB, then it will turn into a row of NA values in the
# taxonomixr::getTaxonomy
taxonomizr_get_taxonomy <- function(
  taxonomy_ids,
  taxonomy_db_path,
  desired_taxa
) {
  result <- withCallingHandlers(
    taxonomizr::getTaxonomy(
      ids = taxonomy_ids,
      sqlFile = taxonomy_db_path,
      desiredTaxa = desired_taxa
    ),
    warning = function(w) {
      # This is the warning that occurs if some of the taxonomy IDs are not
      # present in the DB. But this can occur naturally in rCRUX, so we don't
      # want a scary warning, rather we want a logged message.
      if (grepl("NAs introduced by coercion", w$message, fixed = TRUE)) {
        # TODO: I want to log this with the logger, not with R's warning
        # mechanism.
        #
        # Supress the warning message
        invokeRestart("muffleWarning")
      }
      # Other warnings will propagate normally
    }
  )

  checkmate::assert_matrix(result, nrows = length(taxonomy_ids))

  result
}
