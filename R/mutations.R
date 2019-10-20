#' @export
readVCFMutations <- function(directory, pattern = "*.vcf"){
  parseVCF <- function(x){
    vcf <- vcfR::read.vcfR(x, verbose = FALSE)
    extracted_fields <- vcfR::vcfR2tidy(vcf, single_frame = TRUE, verbose = FALSE)$dat
    extracted_fields <- extracted_fields[,c("CHROM", "POS", "Indiv", "gt_AF")]
    extracted_fields <- extracted_fields[order(extracted_fields$Indiv),]
    return(extracted_fields)
  }
  fileList <- dir(directory, pattern, full.names = TRUE)
  segmentList <- lapply(fileList, parseVCF)
  segmentList <- segmentList[!unlist(lapply(segmentList, function(x){any(is.na(x))}))]
  segmentTable <- data.table::rbindlist(segmentList)
  colnames(segmentTable) <- c("Chr", "Pos", "SampleID", "AF")
  segmentTable$AF <- as.numeric(segmentTable$AF)
  segmentTable <- na.omit(segmentTable)
  return(segmentTable)
}

#' @export
getScoreMutations <- function(pair, segmentTable, populationMutations, nAdditionalSamples = 0){
  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  sample1_granges <- makeGRangesFromDataFrame(sample1, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
  sample2_granges <- makeGRangesFromDataFrame(sample2, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
  # sample1_granges$AF <- as.numeric(sample1_granges$AF)
  # sample2_granges$AF <- as.numeric(sample2_granges$AF)
  sample1_granges$AF <- sample1_granges$AF / max(sample1_granges$AF)
  sample2_granges$AF <- sample2_granges$AF / max(sample2_granges$AF)

  overlaps <- suppressWarnings(findOverlaps(sample1_granges, sample2_granges))
  hits_sample1 <- sample1_granges[queryHits(overlaps)]
  hits_sample2 <- sample2_granges[subjectHits(overlaps)]
  # hits_sample1$AF <- as.numeric(hits_sample1$AF)
  # hits_sample2$AF <- as.numeric(hits_sample2$AF)
  # hits_sample1$AF <- hits_sample1$AF / max(sample1_granges$AF)
  # hits_sample2$AF <- hits_sample2$AF / max(sample2_granges$AF)

  # score <- sum(hits_sample1$AF, hits_sample2$AF) / (0.5 * sum(length(sample1_granges), length(sample2_granges)))

  nSamples <- length(unique(segmentTable$SampleID)) + nAdditionalSamples



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
      (sample1_granges$AF / sqrt(
        countOverlaps(sample1_granges, populationMutations) / nSamples
        )
       ),
      (sample2_granges$AF / sqrt(
        countOverlaps(sample2_granges, populationMutations) / nSamples
        )
       )
      )
    )
  )


  return(score)
}

# getScoreMutations(p, t)

#' @export
getScoresMutations <- function(pairs, segmentTable, additionalMutations = NULL, nAdditionalSamples = 0, reference = NULL, excludeChromosomes = "chrY"){
  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(segmentTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }
  pair_scores <- apply(pairs, 1, function(x){getScoreMutations(as.character(x), segmentTable, populationMutations, nAdditionalSamples)})

  if(is.null(reference)){warning("No reference supplied, p-values not calculated", immediate. = TRUE)}
  pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
  results <- cbind.data.frame(pairs, pair_scores, pair_ps)

  return(results)
}

#' @export
makeReferenceMixingPairsMutations <- function(segmentTable, pairs, nperm = 10, additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y"){

  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(segmentTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }
  reference <- numeric()
  for(i in 1:nperm){
    message("Constructing reference: Iteration #", i)

    random_pairs <- as.data.table(cbind(sample(pairs[[1]]), sample(pairs[[2]])))
    random_pairs <- random_pairs[!apply(random_pairs, 1, function(y){any(apply(pairs, 1, function(x){all(x == y)}))})]
#    print(random_pairs)
    pair_scores <- apply(random_pairs, 1, function(x){getScoreMutations(as.character(x), segmentTable, populationMutations, nAdditionalSamples)})
    reference <- c(reference, pair_scores)


  }
  return(reference)
}

collatePopulationMutations <- function(segmentTable){

  populationMutations <- makeGRangesFromDataFrame(segmentTable[,c("Chr", "Pos")], start.field = "Pos", end.field = "Pos")
  return(populationMutations)
}


