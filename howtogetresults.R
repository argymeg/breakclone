tab <- readAlleleSpecific(directory = "~/Documents/ascat_collate_sloane_finals", pattern = "*Segments_AbsCN_alleleSpecific_profile.txt")
p <- inferPairs(tab)
# ref <- makeReferenceMixingPairs(tab, p, 5)
ref <- readRDS("reference_ASCAT_allrealpairs_run1_200iters.rds")

results <- getScores(p, tab, ref)
plotScores(ref, results)



vcftab <- readVCF("vcf_sample/")
vcfpairs <- data.table::fread("nki_vcf_pairs.txt", header = FALSE)
vcftab$SampleID <- substr(vcftab$SampleID, 12, 16)
# vcfref <- makeReferenceMixingPairs(vcftab, vcfpairs, 200, cnType = "VCF", excludeChromosomes = "chrY")
vcfref <- readRDS("~/Documents/breaklone/chosen_references/reference_VCF_mixpairs_run1_200iters.rds")
resultsvcf <- getScores(vcfpairs, vcftab, vcfref, cnType = "VCF")
plotScores(vcfref, resultsvcf)

tabtaps <- readAlleleSpecific(c("/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/LCIS/", "/Users/argymeg/Documents/lcis-clonality/ALL_TAPS_OUTPUT/100probes/INV/"), pattern = "_segmentCN.txt", nmajor.field = NULL, ntotal.field = "Cn", nminor.field = "mCn", chr.field = "Chromosome", nprobes.field = "probes", sample.field = NULL)
ptaps <- inferPairs(tabtaps)
# reftaps <- makeReferenceMixingPairs(tabtaps, ptaps, 200, excludeChromosomes = "chrY")
reftaps <- readRDS("~/Documents/breaklone/chosen_references/reference_TAPS_100_run1_200iters.rds")
restaps <- getScores(ptaps, tabtaps, excludeChromosomes = "chrY", reference = reftaps)
plotScores(reftaps, restaps)



t <- readVCFMutations("mutations/Clonality_VCFs/")
t <- t[grep("N", t$SampleID, invert = TRUE),]
t$SampleID <- sub("P", "_P", t$SampleID)
t$SampleID <- sub("IR", "_IR", t$SampleID)
p <- inferPairs(t)
p <- p[,2:1]

# library(TCGAretriever)
# library(data.table)
# library(GenomicRanges)
#
# setwd ("/Users/argymeg/Documents/breaklone/mutations")
# BEDFILE <- fread("~/Documents/wgs-mut-clonality/Sloane_Covered_edited.bed")
# GENELIST<-unique(BEDFILE$V4)
# GENELIST<-grep(",", GENELIST, invert = TRUE, value = TRUE)
#
# allmuts <- get_ext_mutation(glist = GENELIST, case_id = "brca_tcga_pan_can_atlas_2018", "brca_tcga_pan_can_atlas_2018_mutations")
# allmuts <- allmuts[allmuts$chr != "NA",]
# allmuts$chr <- paste0("chr", allmuts$chr)
# allmuts$chr <- sub("23", "X", allmuts$chr)
# allmuts <- allmuts[allmuts$start_position == allmuts$end_position,]
# allmuts <- allmuts[allmuts$mutation_type %in% c("Missense_Mutation","Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region", "Splice_Site", "Translation_Start_Site"),]
#
# mutation_freqs <- makeGRangesFromDataFrame(allmuts, start.field = "start_position", end.field = "end_position", seqnames.field = "chr")
#
# allvcfs <- dir(path = "collect_all_vcfs", pattern = "*.vcf", full.names = TRUE, recursive = FALSE)
# for(vcffile in allvcfs){
#   try({
#     avcf <- read.table(vcffile, stringsAsFactors = FALSE)[,1:2]
#     #    avcf$V1 <- sub("chr", "", avcf$V1)
#     avcf <- makeGRangesFromDataFrame(avcf, seqnames.field = "V1", start.field = "V2", end.field = "V2")
#     mutation_freqs <- c(mutation_freqs, avcf)
#   })
# }

mutation_freqs <- readRDS("~/Documents/breaklone/chosen_references/additionalMutations.rds")
# reftcga <- makeReferenceMixingPairsMutations(t, p, 200, additionalMutations = mutation_freqs, nAdditionalSamples = 1161)
reftcga <- readRDS("~/Documents/breaklone/chosen_references/reference_mutations_mixpairs_run1_200iters.rds")
restcga <- getScoresMutations(pairs = p, segmentTable = t, additionalMutations = mutation_freqs, reference = reftcga, nAdditionalSamples = 1161)
plotScores(reftcga, restcga)

