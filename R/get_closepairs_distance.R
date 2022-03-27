#' get_closepairs_distance
#'
#' This function takes a distance matrix and finds all pairs of
#' samples that are within a specified distance cutoff. The pairs are
#' returned as an nx2 matrix, where n is the number of close pairs, and
#' the elements of the matrix are the indices of the samples in that close
#' pair in the distance matrix.
#'
#'@param distance_matrix A distance matrix
#'@param distance_cutoff An inclusive distance cutoff.
#'@param show_hist Boolean indicator for whether or not to plot a histogram of the distances and cutoff. Uses ggplot defaults for bins.
#'@param verbose TRUE to print distance cutoff, FALSE otherwise
#'@export
#'@examples
#'get_closepairs_distance(distance_matrix, 0.0005)

get_closepairs_distance <- function(distance_matrix, distance_cutoff, show_hist, verbose){
  if(verbose){print(paste0('Distance cutoff: ', distance_cutoff))}
  if(show_hist){
    distances <- as.numeric(distance_matrix[upper.tri(distance_matrix)])
    d <- data.frame(distances)
    p <- ggplot(d, aes(x = distances)) + geom_histogram(na.rm = T, bins = 30) +
      theme_bw() + xlab('Pairwise Distances') + ylab('Number of Observations') +
      geom_vline(xintercept = distance_cutoff, lty = 2, color = 'blue', size = 1)
    show(p)
  }
  close_pairs <- arrayInd(which(distance_matrix <= distance_cutoff), dim(distance_matrix)) # This includes self-comparisons and
  # each pair flipped (and therefore double-counted)
  close_pairs <- close_pairs[close_pairs[,1] < close_pairs[,2],] # Removes self-comparisons and one double-count for each pair
  return(close_pairs)
}


