library(GenomicRanges)
library(data.table)

tab <- readAlleleSpecific(c("~/Documents/clonality_newer/salpies_clonality/Primaries/", "~/Documents/clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)

reftab <- readAlleleSpecific("~/Documents/clonality_newer/salpies_clonality/Controls/")
ref <- makeReference(reftab, 5)
ref <- makeReferenceMixingPairs(tab, p, 5)


results <- getScores(p, tab, ref)
# write.csv(results, "ir_still_still_early_results_justthiscohort.csv")
plotScores(ref, results)

vcftab <- readVCF("vcf_sample/")
vcfref <- makeReference(vcftab, 2, "VCF")
randomise <- sample(unique(vcftab$SampleID))
random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
resultsvcf <- getScores(random_pairs, vcftab, vcfref, cnType = "VCF")

tabtaps <- readAlleleSpecific(c("/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/LCIS/", "/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/INV/"), pattern = "_segmentCN.txt", nmajor.field = NULL, ntotal.field = "Cn", nminor.field = "mCn", chr.field = "Chromosome", nprobes.field = "probes", sample.field = NULL)
ptaps <- inferPairs(tabtaps)
tabtaps_ref <- tabtaps[grep("INV", tabtaps$SampleID, invert = TRUE),]
reftaps <- makeReference(tabtaps_ref, 5, excludeChromosomes = "chrY")
reftaps <- makeReferenceMixingPairs(tabtaps, ptaps, 5, excludeChromosomes = "chrY")
restaps <- getScores(ptaps, tabtaps, excludeChromosomes = "chrY", reference = reftaps)
plotScores(reftaps, restaps)


