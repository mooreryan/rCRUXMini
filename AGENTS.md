Functions that are internal to a particular R script should be prefixed with a `.`. E.g.,

```R
# This function is used in other R files
foo <- function() {
  .bar()
}

# This function is only used in this R file
.bar <- function() {
  # ...
}
```
