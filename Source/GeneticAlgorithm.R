
# this genetic algorithm does not make use of a larger 'population' to which 
# the offspring are added to. The offspring become the entire population for 
# the next generation. 
# Parents are selected from the resulting population weighted by a
# performance metric.

# single-point crossover
genetic_algorithm <- function(Theta, scores=rep(1L,nrow(Theta))){
  # default 'scores' corresponds to even weighting across all genes 
  # (score=1 for all genomes)
  g <- nrow(Theta)
  stopifnot(g%%2==0)  #num of genomes must be a multiple of 2 
  p <- ncol(Theta)    #num of genes
  
  Theta_new <- crossover_recombination(Theta, scores) #g and p can be calculated within, but is it cheaper to perhaps calculate once and pass it along each generation?
  
}

crossover_recombination <- function(Theta,scores){
  g <- nrow(Theta)    # num genomes
  stopifnot(g%%2==0)  # num of genomes must be a multiple of 2 
  n <- g/2            # num pairs
  p <- ncol(Theta)    # num genes
  crossmat <- create_crossover_matrix(g,n,p,scores)
  
  # should be able to post/premultiply in a way that gives us the resulting matrix we want
  # PRE multiply
  # each row of the matrix we premult by will represent the target row
  # the column marked with 1 will give the row which will be placed in the target row
  
  parents <- crossmat[,1:2] |> t() |> matrix(ncol=1)
  # transpose to ensure the correct order is used to fill the matrix 
  # (byrow has no effect, it always goes down the column first)
  # now we will have a single column with integers in them.
  # each row is in the correct place, it is just to use the integer
  # to tell us where to place the 1 (all else being 0)
  # if can do this, then we will have the parents in the right place
  # now just to crossover...ahhhh
  
  xindex <- crossmat[,3] #cross-over point index
  
  base <- matrix(0L,ncol=g,nrow=g)
  for (i in 1:g)
  {
    base[i,parents[i]] <- 1L
  }
  selected_genomes <- base %*% Theta
  
  transformed_genomes <- selected_genomes
  for(i in 2*(1:n))
  {
    transformed_genomes[(i-1):i,xindex[i]:p] <- selected_genomes[i:(i-1),xindex[i]:p]
    # in theory, this works well
    # the subset-to-be-flipped is replaced with the (row-flipped) subset 
    # from the original
  }
}

# columns 1:2 represent the parents
# column 3 is the crossover point
# g: num of genomes (parents)
# p: num of genes / parameters within each genome
create_crossover_matrix <- function(g,n,p,scores){
  
  return(
    cbind(
      t(replicate(n, sample(g, size=2, prob=scores))), #selects 2 parents
      sample(p, n, replace=TRUE) #selects crossover point
    ) 
  )
  
  #note: will randomly crossover an amount of 'genes' from 1:num(genes)
  #thus the crossover cannot be 'no crossover', but it can be 'ALL crossover', 
  #resulting in offspring identical to the parents in only one case.
}

# to evaluate the scores applied to each parent based on some (arbitrary)
# performance metric
score_genomes <- function(Theta) {
  scores <- future_pmap(Theta,.f=objective,.progress=TRUE)
  return(scores)
}