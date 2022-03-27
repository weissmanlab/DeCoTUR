#' get_discordances
#'
#' This function takes a presence-absence matrix and a set of close pairs
#' and computes the number of discordances among close pairs for each row of the matrix.
#' Used in the process of obtaining p-values. Will optionally write the discordance information to file.
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample. Rownames should be trait names, colnames should be sample names.
#'@param closepairs A set of close pairs obtained from get_closepairs.
#'@param verbose Boolean to indicate if progress bar is used.
#'@export
#'@examples
#'


get_discordances <- function(pa_matrix, closepairs, verbose){
  if(verbose){
    discdat <- pbapply::pbapply(pa_matrix, 1, function(x){get_discordances_row(x, closepairs)})
  } else{
    discdat <- apply(pa_matrix, 1, function(x){get_discordances_row(x, closepairs)})
  }
  ddiscdat <- data.frame(discdat)
  ddiscdat$trait <- rownames(ddiscdat)
  rownames(ddiscdat) <- NULL
  ddiscdat <- ddiscdat[, c(2, 1)]
  return(ddiscdat)
}

