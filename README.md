_Note: in progress!_

# rCRUXMini

<!-- badges: start -->
<!-- badges: end -->

The goal of rCRUXMini is to ...

## Installation

You can install the development version of rCRUXMini from [GitHub](https://github.com/) with:

```r
# install.packages("pak")
pak::pak("mooreryan/rCRUXMini")
```

## Example

This is a basic example which shows you how to solve a common problem:

```r
library(rCRUXMini)
## basic example code
```

## Hacking Guidelines

### Functions in there own files

This might seem a bit tedious, but putting functions in their own files can help the linter (lintr) to find issues such as `no visible binding for global variable 'apple' (lintr object_usage_linter)`.

This can be especially important when you're using futures with multisession plans, as you have to ensure all dependencies (variables and functions) are made available to the functions executed
in a future. This requires careful use of `.env_globals`, or simply passing in the dependencies as arguments to the function. The second option is much more easily done if your linter can alert you to variables that any of your functions/closures actually close over!

### Mapping to Original rCRUX

- `rCRUX::get_seeds_local` -> `rCRUXMini::pipeline` (we call this finding amplicons)
- `rCRUX::blast_seeds` -> `rCRUXMini::pipeline` (we call these pulling amplicons)
  - They separate out this function into a lot of little blast helpers and each of them do some parsing of the results

I need to separate out the functions as well so users can run them seperately if required

### Assertion Errors

Keep in mind that pretty much all functions may throw a `checkmateError` which signals that an assertion has failed. Ideally these should never be caught, so if you want to be super thorough, you could include a clause specifically for that error in you `try_fetch`/`tryCatch` calls, and rethrow those specific conditions. The reason is because `checkmateError` itself is a sublcass of `error`, which regretably means that catching the `error` condition will also catch assertion errors, and this is almost certainly _not_ what you inteded to do.

---

We need to manually exclude non-R files because lintr in Zed will try to lint non-R files.

---

Plausible amplicon coordinates files have `subject_accession_version`, `subject_gi`, `unique_subject_taxonomy`, etc. This is because those are the subjects from the BLAST DBs that you searched against.

---

If the parent R process is killed...you will get ZOMBIE WORKERS!!!
