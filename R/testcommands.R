library(GenomicRanges)
library(data.table)


tab <- readAlleleSpecific(c("../clonality_newer/salpies_clonality/Primaries/", "../clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)


getScore(as.character(p[4,]), tab)


results <- getScores(p, tab)
write.csv(results, "ir_still_still_early_results_justthiscohort.csv")



reftab <- readAlleleSpecific("../clonality_newer/salpies_clonality/Controls/")


reference <- numeric()
for(i in 1:10){
  print(i)
  randomise <- sample(unique(reftab$SampleID))
  random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
  apply(random_pairs, 1, function(x){if(x[1] == x[2]){stop("yes, it's possible: ", x[1])}})
  reference <- c(reference, getScores(random_pairs, reftab, TRUE))
}
