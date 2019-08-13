readAlleleSpecific <- function(directory, pattern = "*_Segments_AbsCN_alleleSpecific.txt"){
  fileList <- dir(directory, pattern, full.names = TRUE)
  segmentList <- lapply(fileList, data.table::fread)
  segmentTable <- data.table::rbindlist(segmentList)
  return(segmentTable)
}

readVCF <- function(directory, pattern = "*.vcf"){
  parseVCF <- function(x){
    vcf <- vcfR::read.vcfR(x, verbose = FALSE)
    extracted_fields <- vcfR::vcfR2tidy(vcf, info_only = TRUE, info_fields = c("SVTYPE", "END", "BINS", "SVLEN"))$fix[,c("CHROM", "POS", "END", "BINS", "SVTYPE", "SVLEN")]
    samplename <- colnames(vcf@gt)[2]
    return(cbind(samplename, extracted_fields))
  }
  fileList <- dir(directory, pattern, full.names = TRUE)
  segmentList <- lapply(fileList, parseVCF)
  segmentList <- segmentList[!unlist(lapply(segmentList, function(x){any(is.na(x))}))]
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

