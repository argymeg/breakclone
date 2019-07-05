readAlleleSpecific <- function(directory, pattern = "*_Segments_AbsCN_alleleSpecific.txt"){
  fileList <- dir(directory, pattern, full.names = TRUE)
  segmentList <- lapply(fileList, data.table::fread)
  segmentTable <- data.table::rbindlist(segmentList)
  return(segmentTable)
}

inferPairs <- function(segmentTable){
  samples <- unique(segmentTable$SampleID)
  patients <- unique(sub("_.+", "", samples))
  pairs <- lapply(patients, grep, samples)
  pairs <- pairs[lengths(pairs) == 2]
  pairs <- as.data.frame(matrix(unlist(pairs), nrow=length(pairs), byrow=T))
  pairs$V1 <- samples[pairs$V1]
  pairs$V2 <- samples[pairs$V2]
  colnames(pairs) <- c("Sample1", "Sample2")
  return(pairs)
}

