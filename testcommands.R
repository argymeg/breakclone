library(GenomicRanges)
library(data.table)

getScore(as.character(p[4,]), tab)

tab <- readAlleleSpecific(c("/mnt/albyn/argy/salpies_clonality/Primaries/", "/mnt/albyn/argy/salpies_clonality/IR/"))
p <- inferPairs(tab)

reftab <- readAlleleSpecific("/mnt/albyn/argy/salpies_clonality/Controls/")
ref <- makeReference(reftab)

results <- getScores(p, tab)
write.csv(results, "ir_still_still_early_results_justthiscohort.csv")
