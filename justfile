test:
    Rscript --vanilla -e 'devtools::test()'

check:
    Rscript --vanilla -e 'roxygen2::roxygenize(clean = TRUE); devtools::check()'

install: check
    Rscript --vanilla -e 'devtools::install(build_vignettes = TRUE)'
