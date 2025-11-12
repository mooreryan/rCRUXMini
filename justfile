default: test

clean_dll: clear
    Rscript --vanilla -e 'pkgbuild::clean_dll()'

compile_dll: clear
    Rscript --vanilla -e 'pkgbuild::compile_dll()'

clear:
    clear; clear; clear

test: clear roxygenize
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

# clang_tidy:
#     clear;clear;clear; clang-tidy \
#     --header-filter='.*' \
#      --exclude-header-filter='^R.*' \
#     ./src/parse_primer_blast.cpp -- \
#     -I/Library/Frameworks/R.framework/Resources/include \
#     -I/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include \
#     -DNDEBUG -fPIC -Wall -O2

clang_tidy2:
    clear;clear;clear; clang-tidy ./src/parse_primer_blast.cpp -- \
    -DNDEBUG -fPIC -Wall -O2


yo:
    # pkgbuild::clean_dll(); pkgbuild::compile_dll(); devtools::load_all()
    # Sys.time(); unlink("/Users/ryan/Desktop/hiii.tsv"); parse_primer_blast(primer_blast_tsv="/Users/ryan/Projects/uw/rCRUX_world/rCRUXMini/_big_nt_test_data/primer_blast.10000000.tsv", output_tsv="/Users/ryan/Desktop/hiii.tsv", maximum_mismatches=4, minimum_length=150, maximum_length=450, num_threads=8); Sys.time()
    echo hi
