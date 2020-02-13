getScoreMutations <- function(segmentTable, pair, populationMutations, nAdditionalSamples = 0){
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

#' @export
calculateRelatednessMutations <- function(pairs, segmentTable, additionalMutations = NULL, nAdditionalSamples = 0, reference = NULL, excludeChromosomes = "chrY"){
  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(segmentTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }
  pair_scores <- apply(pairs, 1, function(x){getScoreMutations(segmentTable, as.character(x), populationMutations, nAdditionalSamples)})

  if(is.null(reference)){warning("No reference supplied, p-values not calculated", immediate. = TRUE)}
  pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
  results <- cbind.data.frame(pairs, pair_scores, pair_ps)

  return(results)
}

#' #' @export
#' makeReferenceMixingPairsMutations <- function(segmentTable, pairs, nperm = 10, additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y"){
#'
#'   segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
#'   populationMutations <- collatePopulationMutations(segmentTable)
#'   if(!is.null(additionalMutations)){
#'     populationMutations <- c(populationMutations, additionalMutations)
#'   }
#'   reference <- numeric()
#'   for(i in 1:nperm){
#'     message("Constructing reference: Iteration #", i)
#'
#'     random_pairs <- as.data.table(cbind(sample(pairs[[1]]), sample(pairs[[2]])))
#'     random_pairs <- random_pairs[!apply(random_pairs, 1, function(y){any(apply(pairs, 1, function(x){all(x == y)}))})]
#' #    print(random_pairs)
#'     pair_scores <- apply(random_pairs, 1, function(x){getScoreMutations(segmentTable, as.character(x), populationMutations, nAdditionalSamples)})
#'     reference <- c(reference, pair_scores)
#'
#'
#'   }
#'   return(reference)
#' }

#' @export
makeReferenceAllPairsMutations <- function(segmentTable, pairs, patients = NULL, delimiter = "_", additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y"){


  if(is.null(patients)){
    p1 <- sapply(strsplit(pairs$Sample1, delimiter), "[", 1)
    p2 <- sapply(strsplit(pairs$Sample2, delimiter), "[", 1)
    if(all(p1 == p2)){
      patients <- p1
    } else {
      stop("Autodetecting patient IDs failed!")
    }
  }
  patients <- rbind(as.data.table(cbind(patients, pairs$Sample1)), as.data.table(cbind(patients, pairs$Sample2)))
  colnames(patients) <- c("patient", "sample")
  patients <- unique(patients)
  setkey(patients, "sample")

  refPairs <- expand.grid(list(Sample1 = unique(pairs$Sample1), Sample2 = unique(pairs$Sample2)), stringsAsFactors = FALSE)
  refPairs <- refPairs[patients[refPairs$Sample1]$patient != patients[refPairs$Sample2]$patient,]

  message("Making reference based on ", nrow(refPairs), " possible pairs, this might take a while")

  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(segmentTable)
  if(!is.null(additionalMutations)){
    populationMutations <- c(populationMutations, additionalMutations)
  }

  reference <- apply(refPairs, 1, function(x){getScoreMutations(segmentTable, as.character(x), populationMutations, nAdditionalSamples)})


  return(reference)
}


collatePopulationMutations <- function(segmentTable){

  populationMutations <- makeGRangesFromDataFrame(segmentTable[,c("Chr", "Pos")], start.field = "Pos", end.field = "Pos")
  return(populationMutations)
}


