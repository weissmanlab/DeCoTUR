#' qpbinom_modified
#'
#' This function is a modified version of the qpbinom function from the package PoissonBinomial.
#' It brute-forces away a numerical issue in R that the original package did not address.
#' Used in the process of obtaining p-values.
#'
#'@param all same as in PoissonBinomial implementation
#'@export
#'@examples
#'

qpbinom_modified <- function (p, probs, wts = NULL, method = "DivideFFT", lower.tail = TRUE,
          log.p = FALSE)
{
  if (!log.p) {
    if (is.null(p) || any(is.na(p) | p < 0 | p > 1))
      stop("'p' must contain real numbers between 0 and 1!")
  }
  else {
    if (is.null(p) || any(is.na(p) | p > 0))
      stop("'p' must contain real numbers between -Inf and 0!")
  }
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic",
                                "Recursive", "Mean", "GeoMean", "GeoMeanCounter", "Poisson",
                                "Normal", "RefinedNormal"))
  cdf <- ppbinom(NULL, probs, wts, method, lower.tail)
  size <- length(probs)
  if (log.p)
    p <- exp(p)
  n0 <- sum(probs == 0)
  n1 <- sum(probs == 1)
  hi <- size - n0
  range<- n1:hi
  cdf[which(cdf > 0.999)] <- 1
  if (lower.tail)
    Q <- stepfun(cdf[range + 1], c(range, hi), right = TRUE)
  else Q <- stepfun(rev(cdf[range + 1]), c(hi, rev(range)),
                    right = TRUE)
  res <- Q(p)
  res[p == lower.tail] <- hi
  res[p == !lower.tail] <- n1
  return(res)
}
