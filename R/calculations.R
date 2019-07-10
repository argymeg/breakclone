#pair must be a character vector!
getScore <- function(pair, segmentTable, breakpoint_scores){
  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  sample1$nTotal <- sample1$nMajor + sample1$nMinor
  sample2$nTotal <- sample2$nMajor + sample2$nMinor

  sample1$status <- ifelse(sample1$nTotal < 2, "loss", ifelse(sample1$nTotal == 2, "norm", "gain"))
  sample2$status <- ifelse(sample2$nTotal < 2, "loss", ifelse(sample2$nTotal == 2, "norm", "gain"))

  sample1_lists <- c(list(), list(sample1["loss", on = "status"]), list(sample1["norm", on = "status"]), list(sample1["gain", on = "status"]))
  sample2_lists <- c(list(), list(sample2["loss", on = "status"]), list(sample2["norm", on = "status"]), list(sample2["gain", on = "status"]))

  # sample1_lists <- c(list(), list(sample1["loss", on = "status"]), list(sample1["gain", on = "status"]))
  # sample2_lists <- c(list(), list(sample2["loss", on = "status"]), list(sample2["gain", on = "status"]))
  #

  sample1_granges <- lapply(sample1_lists, makeGRangesFromDataFrame)
  sample2_granges <- lapply(sample2_lists, makeGRangesFromDataFrame)

  hits_start <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 10000)))]}, sample1_granges, sample2_granges)
  hits_end <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 10000)))]}, sample1_granges, sample2_granges)

  nonhits_start <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 10000)))]}, sample1_granges, sample2_granges)
  hits_end <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 10000)))]}, sample1_granges, sample2_granges)

  merged_lists <- mapply(function(x, y){list(merge(x, y, by = c("Chr", "Start"), all = TRUE))}, sample1_lists, sample2_lists)
  binary_flags <- lapply(merged_lists, function(x){apply(x, 1, function(y){!any(is.na(y))})})

  concordant_actual_list <- mapply(function(x,y){list(x[y,])}, merged_lists, binary_flags)
  discordant_actual_list <- mapply(function(x,y){list(x[!y,])}, merged_lists, binary_flags)

  #all.x here and na.rm one down are here for dev purposes - missing scores shouldn't be allowed!
  concordant_actual_list_with_scores <- lapply(concordant_actual_list, merge, breakpoint_scores, all.x = TRUE)

  #nconcordant_adj <- sum(unlist(lapply(concordant_actual_list_with_scores, function(x){sum(x$Score, na.rm = TRUE)})))
  nconcordant_adj <- sum(unlist(lapply(concordant_actual_list_with_scores, function(x){sum(x$score)})))
  ndiscordant <- sum(unlist(lapply(discordant_actual_list, nrow)))

  # concordant_list <- lapply(binary_flags, sum)
  # discordant_list <- lapply(binary_flags, function(x){length(x) - sum(x)})
  #
  # nconcordant <- unlist(concordant_list)
  # ndiscordant <- unlist(discordant_list)

  # lapply(ssadistsss_lists, function(x){apply(x, 1, function(y){!any(is.na(y))})})
  # lapply(ssadistsss_lists, function(x){table(apply(x, 1, function(y){any(is.na(y))}))})
  #
  # lapply(sample1_lists, function(s1, s2){lapply(s1$Start, breakpointScore, s2$Start)}, sample2_lists)
  # unlist(mapply(function(s1, s2){lapply(s1$Start, breakpointScore, s2$Start)}, sample1_lists, sample2_lists))


  score <- nconcordant_adj/(nconcordant_adj + 0.5 * ndiscordant)
  return(score)

  # scores <- unlist(lapply(sample1$Start, breakpointScore, sample2$Start))
  # return(scores)
}

# breakpointScore <- function(breakpoint, sample2){
#   if(breakpoint %in% sample2){
#     return(1)
#   } else {
#     return(0)
#   }
# }

calculateBreakpointScores <- function(segmentTable){
  all_breakpoints <- rbind(segmentTable[,c("Chr", "Start")], segmentTable[,c("Chr", "End")], use.names = FALSE)
  breakpoint_scores <- as.data.table(dplyr::count(all_breakpoints, Chr, Start))
#  breakpoint_scores$penalty <- breakpoint_scores$n / length(unique(segmentTable$SampleID))
  breakpoint_scores$score <- 1 - breakpoint_scores$n / length(unique(segmentTable$SampleID))
  return(breakpoint_scores)
}

