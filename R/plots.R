#' @import ggplot2

#' @export
plotScores <- function(reference, results){
  names(reference) <- rep("Reference", length(reference))
  pair_scores <- results$pair_scores
  names(pair_scores) <- rep("Pairs", length(pair_scores))
  toplot <- as.data.table(c(pair_scores, reference), keep.rownames = TRUE)
  names(toplot) <- c("Cohort", "Score")
  plot <- ggplot(toplot, aes(x = Score, colour = Cohort)) + geom_density(size = 1) + scale_colour_brewer(palette = "Set1") +
    geom_hline(yintercept = 0, size = 1) + geom_vline(xintercept = 0, size = 1) + geom_vline(xintercept = quantile(reference, probs = 0.95), size = 0.5, linetype = 2)
  return(plot)
}
