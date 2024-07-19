#' This function modifies a list of simulated datasets by removing the oldest species
#' from random lineages.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages.
#' @param percentage_spec_to_remove Percentage of species to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
 
#' Load necessary packages
library(ape)
library(DAISIE)
library(DDD)



remove_oldest_spec_from_random_lineages <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G) {
  data_sims_G7 <- list()  # Initialize an empty list to store the modified datasets
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees

  # Assuming 'i' is defined somewhere outside this function, as it is used but not defined here
  one_island <- simulated_dataset_to_change_data[[i]]  # Extract the dataset for a specific island
  one_island_tree <- simulated_tree_to_change_data[[i]]  # Extract the corresponding tree data for the same island
  
  # Calculate the number of species to remove based on the given percentage
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
  
  one_island1 <- one_island  # Create a copy of the island dataset to modify
  species_already_removed <- 0  # Initialize the counter for the number of species removed

  # Main loop to continue until the required number of species are removed
  while (species_already_removed < spec_to_remove) {
    # Loop over each lineage in the island dataset, starting from the second element
    for (j in 2:length(one_island1)) {
      # Condition to handle specific cases where the branching times are exactly 2
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max) {
        one_island1[[j]]$branching_times = c(island_age, col_max)  # Update branching times
        species_already_removed <- species_already_removed + 2  # Increment the removed species counter by 2
        
        # Update the status (stac) based on the original stac value
        if (one_island[[j]]$stac == 2) {
          one_island1[[j]]$stac <- 5
        } else if (one_island[[j]]$stac == 3) {
          one_island1[[j]]$stac <- 7
        } else if (one_island[[j]]$stac == 4) {
          one_island1[[j]]$stac <- 1
        }
      }

      # Condition to handle cases where branching times are greater than or equal to 3
      else if (length(one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove) {
        trees <- one_island_tree[[j]]  # Extract the tree for the current lineage
        branching_times <- one_island1[[j]]$branching_times[-c(1,2)]  # Remove the first two branching times
        
        if (length(branching_times) > 1) {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          S <- length(branching_times) + 1
          
          # Calculate the tip lengths and find the branch with the maximum length to drop
          tip_lengths <- setNames(trees$edge.length[sapply(1:S, function(x,y) which(y == x), y = trees$edge[,2])], trees$tip.label)
          branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
          
          if (length(branch_to_drop_max) > 1) {
            branch_to_drop_max <- branch_to_drop_max[1]  # In case of ties, take the first one
          }
          
          # Drop the tip with the maximum length and update the tree
          tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
          tip_lengths <- tip_lengths[-(which(tip_lengths == max(tip_lengths)))[1]]
          trees <- tree_reduced
          
          # Recalculate the branching times and update the dataset
          bt <- branching.times(tree_reduced)
          brts <- sort(as.vector(bt), decreasing = TRUE)
          New_branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
          one_island1[[j]]$branching_times <- New_branching_times
          branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
          
          # Update the removed species count and missing species count
          species_already_removed <- species_already_removed + ((length_old_branching_times) - length(New_branching_times))
          one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + ((length_old_branching_times) - length(New_branching_times))
          one_island_tree[[j]] <- trees
        } else {
          # If there is only one branching time, update the branching times and counters accordingly
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          New_branching_times <- one_island1[[j]]$branching_times[c(1,2)]
          one_island1[[j]]$branching_times <- New_branching_times
          branching_times <- one_island1[[j]]$branching_times[c(1,2)]
          species_already_removed <- species_already_removed + 1
          one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + 1
          one_island_tree[[j]] <- trees
        }
      }

      # Additional condition to handle cases where branching times need to be adjusted
      else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max) {
        if (one_island[[j]]$stac == 2) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        } else if (one_island[[j]]$stac == 3) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 7
        }
      }

      # Final condition to handle cases where the branching times are already at the maximum value
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max) {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
    }

    # After the loop, assign the modified dataset to the output list
    data_sims_G7 <- one_island1
  }

  return (data_sims_G7)
}
