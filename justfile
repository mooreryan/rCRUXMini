default: test

test:
    Rscript --vanilla -e 'devtools::test()'

test_file basename:
    Rscript --vanilla -e 'devtools::test_active_file("tests/testthat/test-{{ basename }}.R")'

check:
    Rscript --vanilla -e 'roxygen2::roxygenize(clean = TRUE); devtools::check()'

install: check
    Rscript --vanilla -e 'devtools::install(build_vignettes = TRUE)'

review_snaps:
    Rscript --vanilla -e 'testthat::snapshot_review()'

accept_snaps:
    Rscript --vanilla -e 'testthat::snapshot_accept()'
