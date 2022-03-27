#' get_scores_pa_closepairs
#'
#' This function takes a presence-absence matrix and a set of close pairs and
#' computes the close-pair score for each pair of objects in the presence-absence matrix (a pair of genes, a pair of snps, etc.).
#'
#'@param pa_matrix A presence-absence matrix for the trait in question (i.e. gene, allele, phenotype, etc.). Rows are trait, columns are sample.
#'@param close_pairs A matrix of close pairs. Output from get_close_pairs_(distance, fraction, auto).
#'@param classweights A vector of weights for each close pair. Output from get_classweights.
#'@param blocksize The number of traits to consider at once in the computation. Changing this number may help speed up the computation. Try 100 first.
#'@param verbose TRUE to print progress, FALSE to not print progress.
#'@param version 'speed' or 'memory' optimized version. 'speed' recommended; if you run into memory issues then decrease the blocksize
#'@export
#'@examples
#'get_scores_pa_closepairs(pa_matrix, close_pairs, classweights, 100, T)

get_scores_pa_closepairs <- function(pa_matrix, close_pairs, classweights, blocksize, verbose, version){
  if(blocksize > dim(pa_matrix)[1]){
    stop("Block size is larger than the number of traits")
  }
  gene_names <- rownames(pa_matrix)
  nrows <- dim(pa_matrix)[1]
  ncols <- dim(pa_matrix)[2]
  numblocks <- (nrows %/% blocksize)
  leftover <- nrows %% blocksize
  if(leftover > 0){
    padding <- blocksize - leftover
    padrows <- matrix(rep(0, padding * ncols), nrow = padding, ncol = ncols)
    rownames(padrows) <- rep('padding', padding)
    colnames(padrows) <- colnames(pa_matrix)
    pa_matrix <- rbind(pa_matrix, padrows)
    numblocks <- numblocks + 1
  }
  scoredat <- c()
  if(verbose){
    pb <- progress_bar$new(
      format = "  Computing scores [:bar] :percent eta: :eta",
      total = numblocks+(leftover > 0),
      clear = FALSE)
    pb$tick(0)
  }
  for(i in 1:numblocks){
    #print(i)
    starti <- (i-1)*blocksize+1
    endi <- i*blocksize
    rangei <- starti:endi
    ind1snp1 <- pa_matrix[rangei, close_pairs[,1]]
    ind2snp1 <- pa_matrix[rangei, close_pairs[,2]]
    for(j in i:numblocks){
      #print(j)
      startj <- (j-1)*blocksize+1
      endj <- j*blocksize
      rangej <- startj:endj
      if(i == j){
        pairindices <- combn(1:blocksize, 2)
      } else{
        pairindices <- t(expand.grid(1:length(rangei), 1:length(rangej)))
      }
      ind1snp2 <- pa_matrix[rangej, close_pairs[,1]]
      ind2snp2 <- pa_matrix[rangej, close_pairs[,2]]
      snpdiffmat1 <- 1*(ind1snp1 == 0 & ind2snp1 == 1) + 2*(ind1snp1 == 1 & ind2snp1 == 0)
      snpdiffmat2 <- 1*(ind1snp2 == 0 & ind2snp2 == 1) + 2*(ind1snp2 == 1 & ind2snp2 == 0)
      if(version == 'speed'){
        test1 <- snpdiffmat1[pairindices[1,],]
        test2 <- snpdiffmat2[pairindices[2,],]
        test12 <- test1 == test2
        test3 <- snpdiffmat1[pairindices[1,],] > 0
        test4 <- snpdiffmat2[pairindices[2,],] > 0
        pos_test <- test12 * test3 * test4
        snpscoremats <- as.numeric(pos_test %*% classweights)
        testn12 <- test1 != test2
        neg_test <- testn12 * test3 * test4
        snpscorematd <- as.numeric(neg_test %*% classweights)
      } else if(version == 'memory'){
        snpscoremats <- pbapply(pairindices, 2, function(x){return((as.numeric(snpdiffmat1[x[1],] == snpdiffmat2[x[2],] & snpdiffmat1[x[1],] > 0 & snpdiffmat2[x[2],] > 0)) %*% classweights)})
        snpscorematd <- pbapply(pairindices, 2, function(x){return((as.numeric(snpdiffmat1[x[1],] != snpdiffmat2[x[2],] & snpdiffmat1[x[1],] > 0 & snpdiffmat2[x[2],] > 0)) %*% classweights)})
      } else {
        stop('Incorrect version. Should be speed or memory')
      }
      snp1 <- rownames(ind1snp1)[pairindices[1,]]
      snp2 <- rownames(ind1snp2)[pairindices[2,]]
      score <- pmax(snpscoremats, snpscorematd)
      scoredat <- rbind(scoredat, data.frame(snp1, snp2, snpscoremats, snpscorematd, score, stringsAsFactors = F))
    }
    if(verbose){pb$tick()}
  }
  names(scoredat) <- c('Trait1', 'Trait2', 'PositiveAssociation', 'NegativeAssociation', 'Score')
  scoredat <- subset(scoredat, Trait1 %in% gene_names & Trait2 %in% gene_names)
  return(scoredat)
}



