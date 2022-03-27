#' get_classweights
#'
#' This function takes a set of close pair classes (output of get_closepair_classes)
#' and computes the weight of each close pair based on the size of the class.
#'@param classes The output of get_closepair_classes. A list of the class of each close pair.
#'@export
#'@examples
#'get_classweights(classes)

get_classweights <- function(classes){
  tc <- table(classes)
  classweights <- sapply(classes, function(x){
    return(1/(tc[which(names(tc) == x)]))
  })
  return(as.numeric(classweights))
}


