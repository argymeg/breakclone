library(GenomicRanges)
library(data.table)


tab <- readAlleleSpecific(c("../clonality_newer/salpies_clonality/Primaries/", "../clonality_newer/salpies_clonality/IR/"))
p <- inferPairs(tab)
#bs <- calculateBreakpointScores(tab)




bigtab <- readAlleleSpecific(c("../clonality_newer/salpies_clonality/Primaries/", "../clonality_newer/salpies_clonality/IR/",
"../clonality_newer/salpies_clonality/DR/","../clonality_newer/salpies_clonality/DR_Contralateral/",
"../clonality_newer/salpies_clonality/IR_Contralateral/", "../clonality_newer/salpies_clonality/Controls/"))


all_breakpoints <- rbind(bigtab[,c("Chr", "Start")], bigtab[,c("Chr", "End")], use.names = FALSE)
abg <- makeGRangesFromDataFrame(all_breakpoints, start.field = "Start", end.field = "Start")


abgs <- makeGRangesFromDataFrame(bigtab[,c("Chr", "Start")], start.field = "Start", end.field = "Start")
abge <- makeGRangesFromDataFrame(bigtab[,c("Chr", "End")], start.field = "End", end.field = "End")


getScore(as.character(p[4,]), tab)

###Don't run!
reference <- numeric()
for(i in 1:10){
  print(i)
  random_pairs <- cbind.data.frame(sample(p$Sample1), sample(p$Sample2))
  reference <- c(reference, apply(random_pairs, 1, function(x){getScore(as.character(x), tab)}))
  print(i)
}

load("reference_10x.RData")


pair_scores <- apply(p, 1, function(x){getScore(as.character(x), tab)})
pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
results <- cbind.data.frame(p, pair_scores, pair_ps)
write.csv(results, "ir_early_results.csv")


