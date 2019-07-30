library(GenomicRanges)
library(data.table)

getScore(as.character(p[4,]), tab)

tab <- readAlleleSpecific(c("../clonality_newer/salpies_clonality/Primaries/", "../clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)

results <- getScores(p, tab)
write.csv(results, "ir_still_still_early_results_justthiscohort.csv")



reftab <- readAlleleSpecific("../clonality_newer/salpies_clonality/Controls/")
ref <- makeReference(reftab)


