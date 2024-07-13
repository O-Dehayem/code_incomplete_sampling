# Load necessary library
library(ape)



# Function to drop the youngest species from a tree
drop_youngest_species <- function(trees, branching_times) {
  
  # Number of species
  S <- length(branching_times) + 1
  
  # Get the tip lengths for each species
  tip_lengths <- setNames(
    trees$edge.length[sapply(1:S, function(x, y) which(y == x), y = trees$edge[, 2])],
    trees$tip.label
  )
  
  # Identify the species with the minimum tip length
  branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
  
  # If there are multiple youngest species, select the first one
  if (length(branch_to_drop_min) > 1) {
    branch_to_drop_min <- branch_to_drop_min[1]
  }
  
  # Drop the identified species from the tree
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
  
  # Get the branching times of the reduced tree
  brts <- branching.times(tree_reduced)
  
  # Sort branching times in decreasing order
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  
  # Return the result as a list
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
}

###############
# Function to drop the oldest species from a tree
drop_oldest_species <- function(trees, branching_times) {
  
  # Number of species
  S <- length(branching_times) + 1
  
  # Get the tip lengths for each species
  tip_lengths <- setNames(
    trees$edge.length[sapply(1:S, function(x, y) which(y == x), y = trees$edge[, 2])],
    trees$tip.label
  )
  
  # Identify the species with the maximum tip length
  branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
  
  # If there are multiple oldest species, select the first one
  if (length(branch_to_drop_max) > 1) {
    branch_to_drop_max <- branch_to_drop_max[1]
  }
  
  # Drop the identified species from the tree
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
  
  # Get the branching times of the reduced tree
  brts <- branching.times(tree_reduced)
  
  # Sort branching times in decreasing order
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  
  # Return the result as a list
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
}
###############
# Function to drop a random species from a tree
drop_random_species <- function(trees, branching_times) {
  
  # Number of species
  S <- length(branching_times) + 1
  
  # Get the tip lengths for each species
  tip_lengths <- setNames(
    trees$edge.length[sapply(1:S, function(x, y) which(y == x), y = trees$edge[, 2])],
    trees$tip.label
  )
  
  # Randomly sample one tip length
  sampled_length <- sample(tip_lengths, 1)
  
  # Get the indices of species with the sampled tip length
  matching_indices <- which(tip_lengths == sampled_length)
  
  # Drop the identified species from the tree
  tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices]))
  
  # Get the branching times of the reduced tree
  brts <- branching.times(tree_reduced)
  
  # Sort branching times in decreasing order
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  
  # Return the result as a list
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
}


