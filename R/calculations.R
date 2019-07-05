#pair must be a character vector!
getScore <- function(pair, segmentTable){
  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  scores <- unlist(lapply(sample1$Start, breakpointScore, sample2$Start))
  return(scores)
}

breakpointScore <- function(breakpoint, sample2){
  if(breakpoint %in% sample2){
    return(1)
  } else {
    return(0)
  }
}
