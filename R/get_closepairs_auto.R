#' get_closepairs_auto
#'
#' This function takes a distance matrix and finds all pairs of
#' samples that are within a distance cutoff determined by an automated algorithm
#' that finds valleys in the histogram of the pairwise distance distribution.
#' The pairs are returned as an nx2 matrix, where n is the number of close pairs, and
#' the elements of the matrix are the indices of the samples in that close
#' pair in the distance matrix. Plotting the histogram with the cutoff is optional.
#' This function was inspired by some code from Lazarus Thurston found at:
#' https://stackoverflow.com/questions/13133297/calculating-peaks-in-histograms-or-density-functions.
#'
#'@param distance_matrix A distance matrix
#'@param which_valley The ordinal valley to be used for the cutoff.
#'@param nbins Number of bins for the histogram
#'@param maxvalleyheight Same as minpeakheight in pracma::findpeaks, which is the function used here. Will probably be need to be adjusted iteratively by the user.
#'@param show_hist Boolean indicator for whether or not to plot a histogram of the distances and cutoff.
#'@param verbose TRUE to report distance cutoff, FALSE to not
#'@export
#'@examples
#'get_closepairs_auto(distance_matrix, 1, 100, T) # This will yield a cutoff for the 1st valley in the pairwise distance distribution,
#'# use 100 bins for the histogram, and shows the plot.

get_closepairs_auto <- function(distance_matrix, which_valley, nbins, maxvalleyheight, show_hist, verbose){
  distances <- as.numeric(distance_matrix[upper.tri(distance_matrix)])
  counts <- hist(distances, breaks = nbins, plot = F)$counts # from Lazarus Thurston
  valleys <- as.numeric(names(counts)[pracma::findpeaks(-counts, minpeakheight = -maxvalleyheight)[,2]])
  if(length(valleys) == 0){
    stop('No valleys detected. Try changing nbins or maxvalleyheight.')
  }
  if(which_valley > length(valleys)){
    stop('There are not enough valleys. Try decreasing which_valley.')
  }
  if(show_plot){
    d <- data.frame(distances)
    p <- ggplot(d, aes(x = distances)) + geom_histogram(bins = nbins, na.rm = T) + geom_vline(xintercept = valleys[-which_valley],lty = 2) +
      theme_bw() + xlab('Pairwise Distances') + ylab('Number of Observations') +
      geom_vline(xintercept = valleys[which_valley], lty = 2, color = 'blue', size = 1)
    show(p)
  }
  distance_cutoff <- valleys[which_valley]
  if(verbose){print(paste0('Distance cutoff: ', distance_cutoff))}
  close_pairs <- arrayInd(which(distance_matrix <= distance_cutoff), dim(distance_matrix)) # This includes self-comparisons and
  # each pair flipped (and therefore double-counted)
  close_pairs <- close_pairs[close_pairs[,1] < close_pairs[,2],] # Removes self-comparisons and one double-count for each pair
  return(close_pairs)
}


