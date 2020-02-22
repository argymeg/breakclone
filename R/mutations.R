getScoreMutations <- function(mutationTable, pair, populationMutations, nAdditionalSamples = 0, scaleAFs){
  sample1 <- mutationTable[mutationTable$SampleID == pair[1],]
  sample2 <- mutationTable[mutationTable$SampleID == pair[2],]

  sample1_granges <- makeGRangesFromDataFrame(sample1, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
  sample2_granges <- makeGRangesFromDataFrame(sample2, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
  # sample1_granges$AF <- as.numeric(sample1_granges$AF)
  # sample2_granges$AF <- as.numeric(sample2_granges$AF)
  if(scaleAFs){
    sample1_granges$AF <- sample1_granges$AF / max(sample1_granges$AF)
    sample2_granges$AF <- sample2_granges$AF / max(sample2_granges$AF)
  }


  overlaps <- suppressWarnings(findOverlaps(sample1_granges, sample2_granges))
  hits_sample1 <- sample1_granges[queryHits(overlaps)]
  hits_sample2 <- sample2_granges[subjectHits(overlaps)]
  if(length(overlaps) > 0){
    nonhits_sample1 <- sample1_granges[-queryHits(overlaps)]
    nonhits_sample2 <- sample2_granges[-subjectHits(overlaps)]
  } else {
    nonhits_sample1 <- sample1_granges
    nonhits_sample2 <- sample2_granges
  }

  # hits_sample1$AF <- as.numeric(hits_sample1$AF)
  # hits_sample2$AF <- as.numeric(hits_sample2$AF)
  # hits_sample1$AF <- hits_sample1$AF / max(sample1_granges$AF)
  # hits_sample2$AF <- hits_sample2$AF / max(sample2_granges$AF)

  # score <- sum(hits_sample1$AF, hits_sample2$AF) / (0.5 * sum(length(sample1_granges), length(sample2_granges)))

  nSamples <- length(unique(mutationTable$SampleID)) + nAdditionalSamples



  score <- sum(
    (hits_sample1$AF + hits_sample2$AF) / sqrt(
      countOverlaps(hits_sample1, populationMutations) / nSamples
      )
    ) / (sum(
      (hits_sample1$AF + hits_sample2$AF) / sqrt(
        countOverlaps(hits_sample1, populationMutations) / nSamples
        )
      ) +
      (0.5 * sum(
      (nonhits_sample1$AF / sqrt(
        countOverlaps(nonhits_sample1, populationMutations) / nSamples
        )
       ),
      (nonhits_sample2$AF / sqrt(
        countOverlaps(nonhits_sample2, populationMutations) / nSamples
        )
       )
      )
    )
  )


  return(score)
}

#' Calculate relatedness scores for paired tumours
#'
#' Calculates the relatedness scores and (optionally) p-values for paired tumours from mutation data
#' @param mutationTable A table of mutations in each sample and their allele frequencies.
#' @param pairs A table of paired samples from the dataset, to test for relatedness.
#' @param additionalMutations A table of mutations to be taken into account when calculating population frequencies. At a minimum, a table of the mutations in the population being studied. This is more informative when tumour type-specific mutations are included from external sources (e.g. TCGA).
#' @param nAdditionalSamples The number of samples used to derive the additional mutations table.
#' @param reference A numeric vector of pair scores comprising the reference distribution, generated from the \code{makeReferenceMutations} function. If omitted, p-value calculation will be skipped.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @return A data frame listing the tumour pairs contained in \code{pairs}, their relatedness scores and p-values for relatedness.

#' @export
calculateRelatednessMutations <- function(mutationTable, pairs, additionalMutations = NULL, nAdditionalSamples = 0, reference = NULL, excludeChromosomes = "chrY", scaleAFs = FALSE){
  mutationTable <- mutationTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(mutationTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }
  pair_scores <- apply(pairs, 1, function(x){getScoreMutations(mutationTable, as.character(x), populationMutations, nAdditionalSamples, scaleAFs)})

  if(is.null(reference)){warning("No reference supplied, p-values not calculated", immediate. = TRUE)}
  pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
  results <- cbind.data.frame(pairs, pair_scores, pair_ps)

  return(results)
}

# #' @export
# makeReferenceMixingPairsMutations <- function(mutationTable, pairs, nperm = 10, additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y"){
#
#   mutationTable <- mutationTable[!excludeChromosomes, on = "Chr"]
#   populationMutations <- collatePopulationMutations(mutationTable)
#   if(!is.null(additionalMutations)){
#     populationMutations <- c(populationMutations, additionalMutations)
#   }
#   reference <- numeric()
#   for(i in 1:nperm){
#     message("Constructing reference: Iteration #", i)
#
#     random_pairs <- as.data.table(cbind(sample(pairs[[1]]), sample(pairs[[2]])))
#     random_pairs <- random_pairs[!apply(random_pairs, 1, function(y){any(apply(pairs, 1, function(x){all(x == y)}))})]
# #    print(random_pairs)
#     pair_scores <- apply(random_pairs, 1, function(x){getScoreMutations(mutationTable, as.character(x), populationMutations, nAdditionalSamples)})
#     reference <- c(reference, pair_scores)
#
#
#   }
#   return(reference)
# }


#' Generate reference distribution from mutation data
#'
#' Generates the reference distribution of concordance scores from unpaired tumours for a given dataset.
#' @param mutationTable A table of mutations in each sample and their allele frequencies.
#' @param pairs A table of paired samples from the dataset. All tumours present in this table will be paired with all tumours from other patients.
#' @param patients A character vector of patient IDs, parallel to the pairs table, used to prevent tumours originating from the same patient from being used in the reference distribution (optional)
#' @param delimiter A character separating patient IDs from tumour-specific identifiers in the sample IDs. Ignored if \code{patients} is provided.
#' @param additionalMutations A table of mutations to be taken into account when calculating population frequencies. At a minimum, a table of the mutations in the population being studied. This is more informative when tumour type-specific mutations are included from external sources (e.g. TCGA).
#' @param nAdditionalSamples The number of samples used to derive the additional mutations table.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @return A numeric vector of pair scores comprising the reference distribution.

#' @export
makeReferenceMutations <- function(mutationTable, pairs, patients = NULL, delimiter = NULL, additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y", scaleAFs = FALSE){
  if(is.null(patients) & is.null(delimiter)){
    patients <- as.character(seq(1, nrow(pairs)))
  } else if(is.null(patients)){
    p1 <- sapply(strsplit(pairs$Sample1, delimiter), "[", 1)
    p2 <- sapply(strsplit(pairs$Sample2, delimiter), "[", 1)
    if(all(p1 == p2)){
      patients <- p1
    } else {
      stop("Autodetecting patient IDs failed!")
    }
  }
  patients <- rbind(as.data.table(cbind(patients, pairs[[1]])), as.data.table(cbind(patients, pairs[[2]])))
  colnames(patients) <- c("patient", "sample")
  patients <- unique(patients)
  setkey(patients, "sample")

  refPairs <- expand.grid(list(Sample1 = unique(pairs[[1]]), Sample2 = unique(pairs[[2]])), stringsAsFactors = FALSE)
  refPairs <- refPairs[patients[refPairs$Sample1]$patient != patients[refPairs$Sample2]$patient,]

  message("Making reference based on ", nrow(refPairs), " possible pairs, this might take a while")

  mutationTable <- mutationTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(mutationTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }

  reference <- apply(refPairs, 1, function(x){getScoreMutations(mutationTable, as.character(x), populationMutations, nAdditionalSamples, scaleAFs)})


  return(reference)
}


collatePopulationMutations <- function(mutationTable){

  populationMutations <- makeGRangesFromDataFrame(mutationTable[,c("Chr", "Pos")], start.field = "Pos", end.field = "Pos")
  return(populationMutations)
}


