library(breaklone)

# tab <- readAlleleSpecific(c("~/Documents/clonality_newer/salpies_clonality/Primaries/", "~/Documents/clonality_newer/salpies_clonality/IR/"))
tab <- readAlleleSpecific(directory = "~/Documents/ascat_collate_sloane_finals", pattern = "*Segments_AbsCN_alleleSpecific_profile.txt")
p <- inferPairs(tab)
# p <- flipPairs(p)

reftab <- readAlleleSpecific("~/Documents/clonality_newer/salpies_clonality/Controls/")
ref <- makeReference(reftab, 5)
ref <- makeReferenceMixingPairs(tab, p, 5)


results <- getScores(p, tab, ref)
# write.csv(results, "ir_still_still_early_results_justthiscohort.csv")
plotScores(ref, results)

vcftab <- readVCF("vcf_sample/")
vcfpairs <- data.table::fread("nki_vcf_pairs.txt", header = FALSE)
vcftab$SampleID <- substr(vcftab$SampleID, 12, 16)
vcfref <- makeReference(vcftab, 2, "VCF", excludeChromosomes = "chrY")
vcfref <- makeReferenceMixingPairs(vcftab, vcfpairs, 5, cnType = "VCF", excludeChromosomes = "chrY")
randomise <- sample(unique(vcftab$SampleID))
random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
resultsvcf <- getScores(vcfpairs, vcftab, vcfref, cnType = "VCF")
plotScores(vcfref, resultsvcf)


tabtaps <- readAlleleSpecific(c("/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/LCIS/", "/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/INV/"), pattern = "_segmentCN.txt", nmajor.field = NULL, ntotal.field = "Cn", nminor.field = "mCn", chr.field = "Chromosome", nprobes.field = "probes", sample.field = NULL)
ptaps <- inferPairs(tabtaps)
tabtaps_ref <- tabtaps[grep("INV", tabtaps$SampleID, invert = TRUE),]
reftaps <- makeReference(tabtaps_ref, 5, excludeChromosomes = "chrY")
reftaps <- makeReferenceMixingPairs(tabtaps, ptaps, 5, excludeChromosomes = "chrY")
restaps <- getScores(ptaps, tabtaps, excludeChromosomes = "chrY", reference = reftaps)
plotScores(reftaps, restaps)


t <- readVCFMutations("mutations/Clonality_VCFs/")
t <- t[grep("N", t$SampleID, invert = TRUE),]
t$SampleID <- sub("P", "_P", t$SampleID)
t$SampleID <- sub("IR", "_IR", t$SampleID)
length(unique(t$SampleID))
p <- inferPairs(t)
p <- p[,2:1]
ref <- makeReferenceMixingPairsMutations(t, p, 20)
res <- getScoresMutations(pairs = p, segmentTable = t, reference = ref)
plotScores(ref, res)

reftcga <- makeReferenceMixingPairsMutations(t, p, 20, additionalMutations = mutation_freqs, nAdditionalSamples = 1161)
restcga <- getScoresMutations(pairs = p, segmentTable = t, additionalMutations = mutation_freqs, reference = reftcga, nAdditionalSamples = 1161)
plotScores(reftcga, restcga)

