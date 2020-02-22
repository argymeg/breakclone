
# breakclone

<!-- badges: start -->
<!-- badges: end -->

Breakclone is a package for assessing clonal relationships of tumour pairs based on copy number or mutation data.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(breakclone)
table <- readAlleleSpecific("/path/to/data")
pairs <- data.table::fread("/path/to/sample/sheet")
reference <- makeReferenceCN(table, pairs)
results <- calculateRelatednessCn(table, pairs, reference)
plotScoresDensity(reference, results)
```

