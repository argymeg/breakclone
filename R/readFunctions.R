#' @import data.table
#'

#' @title Can't believe this
#'
#' @description This either
#'
#' @export
readAlleleSpecific <- function(directory, pattern = "*_Segments_AbsCN_alleleSpecific.txt", sample.field = "SampleID",
                               chr.field = "Chr", start.field = "Start", end.field = "End", nprobes.field = "nProbes",
                               nmajor.field = "nMajor", nminor.field = "nMinor", ntotal.field = NULL){
  fileList <- dir(directory, pattern, full.names = TRUE)
  #segmentList <- lapply(fileList, data.table::fread[,c(sample.field, chr.field, start.field, end.field, nprobes.field, nmajor.field, nminor.field)])
  segmentList <- lapply(fileList, data.table::fread)
  if(is.null(sample.field)){segmentList <- mapply(cbind.data.frame, segmentList, basename(fileList), SIMPLIFY = FALSE)}
  segmentTable <- data.table::rbindlist(segmentList)
  if(is.null(sample.field)){
    sample.field <- "SampleID"
    names(segmentTable)[ncol(segmentTable)] <- "SampleID"
    }

  if(is.null(nmajor.field)){
    nmajor.field = "nMajor"
    segmentTable[,nmajor.field] <- segmentTable[,..ntotal.field] - segmentTable[,..nminor.field]
  }

  cols_needed <- c(sample.field, chr.field, start.field, end.field, nprobes.field, nmajor.field, nminor.field)
  colnames_needed <- c("SampleID", "Chr", "Start", "End", "nProbes", "nMajor", "nMinor")
  segmentTable <- segmentTable[,..cols_needed]
  names(segmentTable) <- colnames_needed


  return(segmentTable)
}

#' @export
readVCF <- function(directory, pattern = "*.vcf"){
  parseVCF <- function(x){
    vcf <- vcfR::read.vcfR(x, verbose = FALSE)
    extracted_fields <- vcfR::vcfR2tidy(vcf, info_only = TRUE, info_fields = c("SVTYPE", "END", "BINS", "SVLEN"))$fix[,c("CHROM", "POS", "END", "BINS", "SVTYPE", "SVLEN")]
    samplename <- colnames(vcf@gt)[2]
    return(cbind(samplename, extracted_fields, stringsAsFactors = FALSE))
  }
  fileList <- dir(directory, pattern, full.names = TRUE)
  segmentList <- lapply(fileList, parseVCF)
  segmentList <- segmentList[!unlist(lapply(segmentList, function(x){any(is.na(x))}))]
  segmentTable <- data.table::rbindlist(segmentList)
  colnames(segmentTable) <- c("SampleID", "Chr", "Start", "End", "Bins", "SVType", "Length")
  return(segmentTable)
}

#' @export
inferPairs <- function(segmentTable, sep = "_"){
  samples <- unique(segmentTable$SampleID)
  patients <- unique(sub(paste0(sep, ".+"), "", samples))
  pairs <- lapply(patients, grep, samples)
  triplets <- pairs[lengths(pairs) == 3]
  pairs <- pairs[lengths(pairs) == 2]
  pairs <- unlist(pairs)
  if(length(triplets) > 0){
    triplets <- as.vector(sapply(triplets, combn, 2))
    pairs <- c(pairs, triplets)
  }
  pairs <- as.data.frame(matrix(pairs, nrow=length(pairs) / 2, byrow=T))
  pairs$V1 <- samples[pairs$V1]
  pairs$V2 <- samples[pairs$V2]
  colnames(pairs) <- c("Sample1", "Sample2")
  return(pairs)
}

#' @export
flipPairs <- function(pairs){
  pairs <- pairs[,c(2,1)]
  return(pairs)
}

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
