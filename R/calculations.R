#pair must be a character vector!
getScore <- function(pair, segmentTable, abg = NULL){

  #This bit should go outside and only happen once!
  all_breakpoints <- rbind(segmentTable[,c("Chr", "Start")], segmentTable[,c("Chr", "End")], use.names = FALSE)
  abg <- makeGRangesFromDataFrame(all_breakpoints, start.field = "Start", end.field = "Start")

  segmentTable <- segmentTable[!"Y", on = "Chr"]

  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  sample1$nTotal <- sample1$nMajor + sample1$nMinor
  sample2$nTotal <- sample2$nMajor + sample2$nMinor

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
  score_from_hits_end <- sum(unlist(lapply(hits_end, function(x){1 - countOverlaps(x, abge, type = "end", maxgap = 5 * 11449) / 208})))

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

# breakpointScore <- function(breakpoint, sample2){
#   if(breakpoint %in% sample2){
#     return(1)
#   } else {
#     return(0)
#   }
# }

# calculateBreakpointScores <- function(segmentTable){
#   all_breakpoints <- rbind(segmentTable[,c("Chr", "Start")], segmentTable[,c("Chr", "End")], use.names = FALSE)
#   breakpoint_scores <- as.data.table(dplyr::count(all_breakpoints, Chr, Start))
# #  breakpoint_scores$penalty <- breakpoint_scores$n / length(unique(segmentTable$SampleID))
#   breakpoint_scores$score <- 1 - breakpoint_scores$n / length(unique(segmentTable$SampleID))
#   return(breakpoint_scores)
# }

