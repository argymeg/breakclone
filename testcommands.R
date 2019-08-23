library(GenomicRanges)
library(data.table)

tab <- readAlleleSpecific(c("~/Documents/clonality_newer/salpies_clonality/Primaries/", "~/Documents/clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)

reftab <- readAlleleSpecific("~/Documents/clonality_newer/salpies_clonality/Controls/")
ref <- makeReference(reftab, 5)

results <- getScores(p, tab)
# write.csv(results, "ir_still_still_early_results_justthiscohort.csv")


vcftab <- readVCF("vcf_sample/")
vcfref <- makeReference(vcftab, 2, "VCF")
randomise <- sample(unique(vcftab$SampleID))
random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
resultsvcf <- getScores(random_pairs, vcftab, vcfref, cnType = "VCF")

tabtaps <- readAlleleSpecific("/Users/argymeg/Downloads/DCIS_Sloane/data/SCINS/materials/tables/data/TAPS", pattern = "_segmentCN.txt", nmajor.field = NULL, ntotal.field = "Cn", nminor.field = "mCn", chr.field = "Chromosome", nprobes.field = "probes", sample.field = NULL)
ptaps <- inferPairs(tabtaps)
tabtaps_ref <- tabtaps[grep("exn", tabtaps$SampleID),]
reftaps <- makeReference(tabtaps_ref, 5)
restaps <- getScores(ptaps, tabtaps, excludeChromosomes = "chrY", reference = reftaps)
