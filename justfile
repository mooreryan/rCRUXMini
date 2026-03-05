default: test

clear:
    clear; clear; clear

test: clear roxygenize install_scripts
    Rscript --vanilla -e 'devtools::test()'

test_file basename: clear
    Rscript --vanilla -e 'devtools::test_active_file("tests/testthat/test-{{ basename }}.R")'

check:
    Rscript --vanilla -e 'roxygen2::roxygenize(clean = TRUE); devtools::check()'

roxygenize:
    Rscript --vanilla -e 'roxygen2::roxygenize(clean = TRUE)'

install: check
    Rscript --vanilla -e 'devtools::install(build_vignettes = TRUE)'

review_snaps:
    Rscript --vanilla -e 'testthat::snapshot_review()'

accept_snaps:
    Rscript --vanilla -e 'testthat::snapshot_accept()'

build_scripts:
    cd scripts && cargo build --release

install_scripts:
    cd scripts && sudo cargo install --path . --root /usr/local --force

instrument_scripts:
    cd scripts && time cargo instruments --release -t time --bin rCRUXMini__ParsePrimerBlast -- --out-file HI --threads=4 --maximum-mismatches=4 --minimum-length=150 --maximum-length=450 ~/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.5000000.tsv
