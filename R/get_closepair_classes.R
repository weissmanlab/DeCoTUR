#' get_closepairs_classes
#'
#' This function takes a set of close pairs and sorts them into bushes (equivalence classes defined by the relation
#' of being a close pair). So if Samples 1 and 2 are close pairs and Samples 2 and 3 are close pairs, then Samples 1, 2, and 3
#' are all in the same close pair class.
#'
#'@param close_pairs A matrix of close pairs; the output of get_closepairs_(distance, fraction, auto)
#'@export
#'@examples
#'get_closepair_classes(close_pairs)

get_closepair_classes <- function(close_pairs){
  relations <- igraph::graph_from_data_frame(close_pairs, directed = F)
  comps <- igraph::components(relations)
  memberships <- comps$membership
  names(memberships) <- names(igraph::V(relations))
  classes <- sapply(close_pairs[,1], function(x){
    classid <- as.numeric(memberships[which(names(memberships) == x)])
    return(classid)
  })
  return(classes)
}


