# This is a test script that makes sure my code is computing the score correctly.

# 5 individuals
n = 5
ng <- 4
# Here's the distance matrix
set.seed(1)
distance_matrix <- matrix(runif(n^2, 0, 1), nrow = 5, ncol = 5)
distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]
diag(distance_matrix) <- 0
distance_cutoff <- 0.6
closepairs <- get_closepairs_distance(distance_matrix, distance_cutoff, T, T)
# There should be 6 close pairs
classes <- get_closepair_classes(closepairs)
# let's get rid of the 1,4 and 2,5
distance_matrix[1,4] <- 0.7
distance_matrix[4,1] <- 0.7
distance_matrix[2,5] <- 0.7
distance_matrix[5,2] <- 0.7
closepairs <- get_closepairs_distance(distance_matrix, distance_cutoff, T, T)
# Don't actually need that for this exercise.
closepairs <- rbind(c(1, 2), c(1, 3), c(2, 3), c(4, 5))
classes <- get_closepair_classes(closepairs)
classweights <- get_classweights(classes)
# let's generate the gene vectors for these guys. 10 genes
generate_gene_vec <- function(ng, seed, thresh){
  set.seed(seed)
  return(as.numeric(runif(ng, 0, 1) < thresh))
}
gene_matrix <- sapply(1:n, function(x){generate_gene_vec(ng, x, 0.5)})
rownames(gene_matrix) <- LETTERS[1:ng]
# OK. The close pairs are 1/2/3 and 4/5
scores <- get_scores_pa_closepairs(gene_matrix, closepairs, classweights, 3, T, 'speed')
# This is correct.

