#' get_discordances_row
#'
#' This function takes a single row of a presence-absence matrix and a set of close pairs
#' and computes the number of discordances among close pairs for that row. Used in
#' the process of obtaining p-values. Helper function for get_discordances.
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample. Rownames should be trait names, colnames should be sample names.
#'@param closepairs A set of close pairs obtained from get_closepairs.
#'@export
#'@examples
get_discordances_row <- function(pa_matrix_row, closepairs){
  name <- rownames(pa_matrix_row)
  pa_vector <- t(as.numeric(pa_matrix_row))
  pa_states <- cbind(pa_vector[closepairs[,1]], pa_vector[closepairs[,2]])
  discordances <- sum(pa_states[,1] != pa_states[,2])
  return(c(name, discordances))
}
