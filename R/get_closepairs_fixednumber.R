#' get_closepairs_fixednumber
#'
#' This function takes a distance matrix and finds all pairs of
#' samples that are within a distance cutoff determined by a specified number of close pairs.
#' If there are too many close pairs that could be below this threshold, then a
#' random sample of these close pairs is returned. User must specify a seed.
#' The pairs are returned as an nx2 matrix, where n is the number of close pairs, and
#' the elements of the matrix are the indices of the samples in that close
#' pair in the distance matrix. Plotting the histogram with the cutoff is optional.
#'
#'@param distance_matrix A distance matrix
#'@param number_close_pairs The number of close pairs to be returned
#'@param seed Random seed for down-sampling if necessary
#'@param show_hist Boolean indicator for whether or not to plot a histogram of the distances and cutoff. Uses ggplot defaults for bins.
#'@param verbose TRUE for returning distance cutoff used, FALSE for not
#'@export
#'@examples
#'get_closepairs_fixednumber(distance_matrix, 5000, 1, T) # This will yield the 5000 closest pairs, with random seed 1.

get_closepairs_fixednumber <- function(distance_matrix, number_close_pairs, seed, show_hist, verbose){
  distances <- as.numeric(distance_matrix[upper.tri(distance_matrix)])
  number_samples <- dim(distance_matrix)[1]
  total_pairs <- (number_samples)*(number_samples - 1)/2
  distance_cutoff <- as.numeric(quantile(distances, number_close_pairs/total_pairs))
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
  if(dim(close_pairs)[1] > number_close_pairs){
    if(verbose){print(paste0('Downsampling ', dim(close_pairs)[1], ' close pairs to ', number_close_pairs, ' close pairs.' ))}
    set.seed(seed)
    samplevec <- sample(1:dim(close_pairs)[1], number_close_pairs,replace = F)
    close_pairs_toreturn <- close_pairs[samplevec,]
    close_pairs <- close_pairs_toreturn
  }
  return(close_pairs)
}


