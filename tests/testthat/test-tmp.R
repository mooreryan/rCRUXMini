# log <- function(msg) {
#   msg <- sprintf("%s -- %s", date(), msg)
#   print(msg)
# }

# describe("thing", {
#   it("does something", {
#     n_max <- 1000000

#     log("reading data")
#     primer_blast_data <- readr::read_tsv(
#       "/Users/ryan/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.tsv",
#       col_types = "ccciiic",
#       n_max = n_max
#     )
#     log("reading data (DONE)")

#     cpp_outfile <- "/Users/ryan/Desktop/yo.tsv"
#     log("C++ parse_primer_blast")
#     parse_primer_blast(
#       primer_blast_tsv = sprintf(
#         "/Users/ryan/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.%d.tsv",
#         n_max
#       ),
#       output_tsv = cpp_outfile,
#       maximum_mismatches = 4,
#       minimum_length = 150,
#       maximum_length = 450
#     )
#     log("C++ parse_primer_blast (DONE)")

#     log("parse_primer_blast_results")
#     primer_blast_results <- parse_primer_blast_results(
#       primer_blast_data,
#       maximum_mismatches = 4
#     )
#     log("parse_primer_blast_results (DONE)")

#     log("find_plausible_amplicon_coordinates")
#     orig <- find_plausible_amplicon_coordinates(
#       primer_blast_results,
#       minimum_length = 150,
#       maximum_length = 450,
#       maximum_mismatches = 4
#     )
#     log("find_plausible_amplicon_coordinates (DONE)")
#   })
# })

# # describe("big data is slow", {
# #   it("parse_primer_blast_results", {
# #     n_max <- 150000

# #     print("---")
# #     log("READING DATA")
# #     primer_blast_data <- readr::read_tsv(
# #       "/Users/ryan/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.tsv",
# #       col_types = "ccciiic",
# #       n_max = n_max
# #       # n_max = 25
# #     )
# #     log("READING DATA (DONE)")

# #     # print("------------------------------------------- ORIGINAL")

# #     # log("PARSING DATA")
# #     # primer_blast_results1 <- parse_primer_blast_results__original(
# #     #   primer_blast_data,
# #     #   maximum_mismatches = 4
# #     # )
# #     # log("PARSING DATA (DONE)")

# #     print("------------------------------------------- DTPLYR")

# #     log("PARSING DATA (B)")
# #     primer_blast_results2 <- parse_primer_blast_results(
# #       primer_blast_data,
# #       maximum_mismatches = 4
# #     )
# #     log("PARSING DATA (DONE)")

# #     # expect_equal(primer_blast_results1, primer_blast_results2)

# #     # ---------------

# #     # log("reading and doing the c++ versoin")
# #     # tmp <- "/Users/ryan/Desktop/tmp_z.tsv"
# #     # yo_parse_primer_blast_results(
# #     #   "/Users/ryan/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.10000000.tsv",
# #     #   tmp,
# #     #   4
# #     # )
# #     # log("reading and doing the c++ versoin (DONE)")
# #     # primer_blast_results3 <- readr::read_tsv(tmp, col_types = "ccciiic")

# #     # expect_equal(
# #     #   primer_blast_results1 |>
# #     #     dplyr::arrange(dplyr::across(dplyr::everything())),
# #     #   primer_blast_results3 |>
# #     #     dplyr::arrange(dplyr::across(dplyr::everything()))
# #     # )

# #     log("find_plausible_amplicon_coordinates")
# #     primer_blast_results2 |>
# #       find_plausible_amplicon_coordinates(
# #         minimum_length = 150,
# #         maximum_length = 450,
# #         maximum_mismatches = 4
# #       )
# #     log("find_plausible_amplicon_coordinates (DONE)")
# #   })
# # })
