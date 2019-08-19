library(GenomicRanges)
library(data.table)

tab <- readAlleleSpecific(c("~/Documents/clonality_newer/salpies_clonality/Primaries/", "~/Documents/clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)

reftab <- readAlleleSpecific("~/Documents/clonality_newer/salpies_clonality/Controls/")
ref <- makeReference(reftab, 2)

results <- getScores(p, tab)
# write.csv(results, "ir_still_still_early_results_justthiscohort.csv")


vcftab <- readVCF("vcf_sample/")
vcfref <- makeReference(vcftab, 2, "VCF")
randomise <- sample(unique(vcftab$SampleID))
random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
resultsvcf <- getScores(random_pairs, vcftab, vcfref, cnType = "VCF")
