#' get_scores
#'
#' This function takes a presence-absence matrix, a distance matrix, a choice of
#' how to determine close pairs, and a list of relevant parameters for that choice,
#' and returns the close pair scores. There is an option for whether or not to
#' downweight bushes. There is also an option to return p-values, in which case three
#' sets of p-values are returned: uncorrected, Bonferroni corrected, and BH corrected.
#' Finally there is an option to run (relatively) speed-optimized or memory-optimized versions of
#' this process.
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample. Rownames should be trait names, colnames should be sample names.
#'@param distance_matrix A distance matrix. Rownames and colnames should be sample names and should be in the same order as the columns in pa_matrix.
#'@param closepair_method Either 'distance', 'fraction', 'fixednumber', or 'auto'.
#'@param closepair_params A list. For 'distance', the distance cutoff and show_hist. For 'fraction', the fraction and show_hist. For 'fixednumber', the number of close pairs, a random seed, and show_hist. For 'auto', the which_valley, nbins, maxvalleyheight, and show_hist. verbose is inherited from the main function call.
#'@param blocksize The number of traits to consider at once in the computation. Changing this number may help speed up the computation. Must be less than or equal to the number of traits.
#'@param downweight TRUE to downweight bushes. Default TRUE.
#'@param withpval TRUE to return p-values. Default TRUE.
#'@param verbose TRUE for progress reports. Default TRUE.
#'@param version 'speed' for speed-optimized, 'memory' for memory optimized. Default 'speed'; if you run into memory problems, try decreasing blocksize first.
#'@export
#'@examples
#'get_scores(pa_matrix, distance_matrix, 'fraction', list(0.1, TRUE), 10, TRUE, TRUE, TRUE)
#'

get_scores <- function(pa_matrix, distance_matrix, closepair_method, closepair_params, blocksize, downweight = TRUE, withpval = TRUE, verbose = TRUE, version = 'speed'){
  if(verbose){print('Starting function.')}
  if(closepair_method == 'distance'){
    if(length(closepair_params) != 2){
      stop('Incorrect number of closepair_params. (Should be 2).')
    }
    distance_cutoff <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_distance(distance_matrix, distance_cutoff, show_hist, verbose)
  } else if(closepair_method == 'fraction'){
    if(length(closepair_params) != 2){
      stop('Incorrect number of closepair_params. (Should be 2).')
    }
    distance_fraction <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_fraction(distance_matrix, distance_fraction, show_hist, verbose)
  } else if(closepair_method == 'auto'){
    if(length(closepair_params) != 4){
      stop('Incorrect number of closepair_params. (Should be 4).')
    }
    which_valley <- closepair_params[[1]]
    nbins <- closepair_params[[2]]
    maxvalleyheight <- closepair_params[[3]]
    show_hist <- closepair_params[[4]]
    close_pairs <- get_closepairs_auto(distance_matrix, which_valley, nbins, maxvalleyheight, show_hist, verbose)
  } else if(closepair_method == 'fixednumber'){
    number_close_pairs <- closepair_params[[1]]
    seed <- closepair_params[[2]]
    show_hist <- closepair_params[[3]]
    close_pairs <- get_closepairs_fixednumber(distance_matrix, number_close_pairs, seed, show_hist, verbose)
  } else{stop('Unidentified close pair method (should be one of: distance, fraction, fixednumber, auto).')}
  if(verbose){
    print(paste0('Obtained ',  dim(close_pairs)[1], ' close pairs.'))
  }
  classes <- get_closepair_classes(close_pairs)
  if(verbose){print('Obtained close pair classes.')}
  if(downweight){
    class_weights <- get_classweights(classes)
  } else{
    class_weights <- rep(1, length(classes))
  }
  if(verbose){print('Obtained close pair class weights.')}
  scores <- get_scores_pa_closepairs(pa_matrix, close_pairs, class_weights, blocksize, verbose, version)
  if(verbose){print('Obtained scores.')}
  if(withpval){
    if(verbose){print('Computing discordance information.')}
    discdat <- get_discordances(pa_matrix, close_pairs, verbose)
    if(verbose){print('Obtained discordance information.')}
    if(downweight){
      unweighted_class_weights <- rep(1, length(classes))
      if(verbose){print('Computing unweighted scores for p-value purposes.')}
      unweighted_scores <- get_scores_pa_closepairs(pa_matrix, close_pairs, unweighted_class_weights, blocksize, verbose, version)
      if(verbose){print('Obtained unweighted scores for p-value purposes.')}
    } else{
      unweighted_scores <- scores
    }
    if(verbose){print('Computing p values')}
    unweighted_scores$sum <- unweighted_scores$PositiveAssociation + unweighted_scores$NegativeAssociation
    unweighted_scores$disc1 <- discdat$disc[match(unweighted_scores$Trait1, discdat$trait)]
    unweighted_scores$disc2 <- discdat$disc[match(unweighted_scores$Trait2, discdat$trait)]
    numclosepairs <- dim(close_pairs)[1]
    fex <- data.frame(unweighted_scores$sum, unweighted_scores$disc1 - unweighted_scores$sum, unweighted_scores$disc2 - unweighted_scores$sum, numclosepairs - unweighted_scores$disc1 - unweighted_scores$disc2 + unweighted_scores$sum)
    if(verbose){
      fres <- pbapply::pbapply(fex,1, function(x) fisher.test(matrix(x,nr=2), alternative = 'greater')$p.value)
    } else{
      fres <- apply(fex,1, function(x) fisher.test(matrix(x,nr=2), alternative = 'greater')$p.value)
    }
    scores$Pval <- fres
    scores$Pval_Bonferroni <- p.adjust(fres, method = 'bonferroni')
    scores$Pval_BH <- p.adjust(fres, method = 'BH')
    if(verbose){print('Obtained p-values')}
  }
  return(scores)
}

