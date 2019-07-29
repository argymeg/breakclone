getScore <- function(pair, segmentTable, abgs = abgs, abge = abge){

  #This bit should go outside and only happen once!
  #all_breakpoints <- rbind(segmentTable[,c("Chr", "Start")], segmentTable[,c("Chr", "End")], use.names = FALSE)
  #abg <- makeGRangesFromDataFrame(all_breakpoints, start.field = "Start", end.field = "Start")

  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  sample1$nTotal <- sample1$nMajor + sample1$nMinor
  sample2$nTotal <- sample2$nMajor + sample2$nMinor

  sample1 <- handlePloidy(sample1)
  sample2 <- handlePloidy(sample2)

  sample1$status <- ifelse(sample1$nTotal < 2, "loss", ifelse(sample1$nTotal == 2 & sample1$nMinor == 1, "norm", ifelse(sample1$nTotal == 2 & sample1$nMinor == 0, "cnloh", ifelse(sample1$nTotal > 4, "amp", "gain"))))
  sample2$status <- ifelse(sample2$nTotal < 2, "loss", ifelse(sample2$nTotal == 2 & sample2$nMinor == 1, "norm", ifelse(sample2$nTotal == 2 & sample2$nMinor == 0, "cnloh", ifelse(sample2$nTotal > 4, "amp", "gain"))))

  sample1 <- sample1[!"norm", on = "status"]
  sample2 <- sample2[!"norm", on = "status"]


  sample1_lists <- c(list(), list(sample1["loss", on = "status"]), list(sample1["cnloh", on = "status"]), list(sample1["gain", on = "status"]), list(sample1["amp", on = "status"]), list(sample1["gain", on = "status"]), list(sample1["loss", on = "status"]))
  sample2_lists <- c(list(), list(sample2["loss", on = "status"]), list(sample2["cnloh", on = "status"]), list(sample2["gain", on = "status"]), list(sample2["amp", on = "status"]), list(sample2["amp", on = "status"]), list(sample2["cnloh", on = "status"]))

  # sample1_lists <- c(list(), list(sample1["loss", on = "status"]), list(sample1["gain", on = "status"]))
  # sample2_lists <- c(list(), list(sample2["loss", on = "status"]), list(sample2["gain", on = "status"]))
  #

  # sample1_granges <- lapply(sample1_lists, makeGRangesFromDataFrame)
  # sample2_granges <- lapply(sample2_lists, makeGRangesFromDataFrame)

  #tryCatch creates an empty GRanges object if the list is empty - would error out otherwise
  sample1_granges <- lapply(sample1_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})
  sample2_granges <- lapply(sample2_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})

  hits_start <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 5 * 11449)))]}, sample1_granges, sample2_granges)
  hits_end <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 5 * 11449)))]}, sample1_granges, sample2_granges)

  score_from_hits_start <- sum(unlist(lapply(hits_start, function(x){1 - countOverlaps(x, abgs, type = "start", maxgap = 5 * 11449) / 208})))
  if(score_from_hits_start < 0){
    warning("Hit next probe! Consider lowering the maxgap")
    score_from_hits_start <- 0
    }
  score_from_hits_end <- sum(unlist(lapply(hits_end, function(x){1 - countOverlaps(x, abge, type = "end", maxgap = 5 * 11449) / 208})))
  if(score_from_hits_end < 0){
    warning("Hit next probe! Consider lowering the maxgap")
    score_from_hits_end <- 0
    }

  # nonhits_start_1 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 11449)))]}, sample1_granges, sample2_granges)
  # nonhits_end_1 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 11449)))]}, sample1_granges, sample2_granges)
  #
  # nonhits_start_2 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 11449)))]}, sample2_granges, sample1_granges)
  # nonhits_end_2 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 11449)))]}, sample2_granges, sample1_granges)

  nconcordant_adj <- score_from_hits_start + score_from_hits_end
  #ndiscordant <- 2 * (nrow(sample1) + nrow(sample2)) - sum(unlist(lapply(c(hits_start, hits_end), length)))
  ndiscordant <- nrow(sample1) + nrow(sample2) - 2 * nconcordant_adj

  score <- nconcordant_adj/(nconcordant_adj + 0.5 * ndiscordant)
  return(score)
}

handlePloidy <- function(sample){
  sample$segLen <- sample$End - sample$Start
  sample_states <- aggregate(segLen ~ nTotal, sample, sum)
  sample_states$adjLen <- as.numeric(sample_states$nTotal) * sample_states$segLen
  sample_ploidy <- sum(sample_states$adjLen) / sum(sample$segLen)

  if(sample_ploidy >= 3.5){
    warning("Sample ", sample$SampleID[1], " looks WGD, inferring ancestral Cn state")
    sample$nTotal <- ifelse(sample$nMinor == (sample$nTotal / 2), 2, sample$nTotal)
    sample$nTotal <- ifelse(sample$nMinor == 0, 1, sample$nTotal)
    sample$nTotal <- ifelse(sample$nTotal >= 3 & sample$nTotal <= 4 & sample$nMinor == 1, 2, sample$nTotal)
    sample$nTotal <- ifelse(sample$nTotal >= 5 & sample$nTotal <= 10, 3, sample$nTotal)
  }
  return(sample)
}

getScores <- function(pairs, segmentTable, isRef = FALSE){
  segmentTable <- segmentTable[!"Y", on = "Chr"]

  filteredTable <- segmentTable
  filteredTable <- filteredTable[!(filteredTable$nMajor == 1 & filteredTable$nMinor == 1),]

  abgs <- makeGRangesFromDataFrame(filteredTable[,c("Chr", "Start")], start.field = "Start", end.field = "Start")
  abge <- makeGRangesFromDataFrame(filteredTable[,c("Chr", "End")], start.field = "End", end.field = "End")

  pair_scores <- apply(pairs, 1, function(x){getScore(as.character(x), segmentTable, abgs, abge)})

  if(isRef){
    results <- pair_scores
  } else {
    pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
    results <- cbind.data.frame(pairs, pair_scores, pair_ps)
  }

  return(results)
}


