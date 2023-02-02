#' filter_traits_by_closepairs
#'
#' This function takes a presence-absence matrix and a set of close pairs and
#' keeps only traits that differ more than once across the set of samples
#' in the close pairs. It also keeps only samples that are present in the
#' closepairs. it also re-indices the close pairs to match the new matrix
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample. Rownames should be trait names, colnames should be sample names.
#'@param closepairs A matrix of close pairs, with each row a pair and each column one of the two sample indices
#'@export
#'@examples


filter_traits_by_closepairs <- function(pa_matrix, closepairs){
  closepairinds <- sort(unique(c(closepairs[,1], closepairs[,2])))
  pa_matrixsub <- pa_matrix[,closepairinds]
  pasum <- rowSums(pa_matrixsub)
  whichpa <- which(pasum > 1 & pasum < length(closepairinds))
  # we need new closepairs which correspond to the new matrix
  cp1 <- match(closepairs[,1], closepairinds)
  cp2 <- match(closepairs[,2], closepairinds)
  newclosepairs <- cbind(cp1, cp2)
  return(list(pa_matrixsub[whichpa,], newclosepairs))
}
