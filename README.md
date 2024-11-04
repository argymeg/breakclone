> [!NOTE]  
> This repository contains the original `breakclone` package described in:
> 
> Lips, E.H., Kumar, T., Megalios, A. *et al*. Genomic analysis defines clonal relationships of ductal carcinoma in situ and recurrent invasive breast cancer. *Nat Genet* **54**, 850–860 (2022). https://doi.org/10.1038/s41588-022-01082-3
> 
> It is no longer updated or maintained. Ongoing development is hosted at: https://github.com/Sawyer-s-Group/breakclone.
# breakclone

<!-- badges: start -->
<!-- badges: end -->

Breakclone is a package for assessing clonal relationships of tumour pairs based on copy number or mutation data.

## Example

This is a basic example of usage with copy number data:

``` r
library(breakclone)
table <- readAlleleSpecific("/path/to/data")
pairs <- data.table::fread("/path/to/sample/sheet")
reference <- makeReferenceCN(table, pairs)
results <- calculateRelatednessCn(table, pairs, reference)
plotScoresDensity(reference, results)
```

This is a very similar example with mutation data:

``` r
library(breakclone)
table <- readVCFMutations("/path/to/data")
pairs <- data.table::fread("/path/to/sample/sheet")
reference <- makeReferenceMutations(table, pairs)
results <- calculateRelatednessMutations(table, pairs, reference)
plotScoresDensity(reference, results)
```
