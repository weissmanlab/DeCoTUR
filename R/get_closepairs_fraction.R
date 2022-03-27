#' get_closepairs_fraction
#'
#' This function takes a distance matrix and finds all pairs of
#' samples that are within a distance cutoff determined by a specified fraction of distances.
#' The pairs are returned as an nx2 matrix, where n is the number of close pairs, and
#' the elements of the matrix are the indices of the samples in that close
#' pair in the distance matrix. Plotting the histogram with the cutoff is optional.
#'
#'@param distance_matrix A distance matrix
#'@param distance_fraction The fraction of distances to be included in the close pairs.
#'@param show_hist Boolean indicator for whether or not to plot a histogram of the distances and cutoff. Uses ggplot defaults for bins.
#'@param verbose TRUE to show distance cutoff, FALSE to not
#'@export
#'@examples
#'get_closepairs_fraction(distance_matrix, 0.1) # This will yield all pairs that are in the smallest 10% for pairwise distance.

get_closepairs_fraction <- function(distance_matrix, distance_fraction, show_hist, verbose){
  distances <- as.numeric(distance_matrix[upper.tri(distance_matrix)])
  distance_cutoff <- as.numeric(quantile(distances, distance_fraction))
  if(verbose){print(paste0('Distance cutoff: ', distance_cutoff))}
  if(show_hist){
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


