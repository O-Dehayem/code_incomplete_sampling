#packages needed
library(ape)



drop_youngest_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
  # Check if branch_to_drop_min has multiple values
  if (length(branch_to_drop_min) > 1) {
    #Drop only one of the youngest species
    branch_to_drop_min <- branch_to_drop_min[1]
  }
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}

###############
drop_oldest_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
  # Check if branch_to_drop_min has multiple values
  if (length(branch_to_drop_max) > 1) {
    # Drop only one of the youngest species
    branch_to_drop_max <- branch_to_drop_max[1]
  }
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}

###############
drop_random_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  sampled_length <- sample(tip_lengths, 1)
  
  # Get the indices of elements equal to the sampled length
  matching_indices <- which(tip_lengths == sampled_length)
  
  # drop the species
  
  
  tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices]))
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}


