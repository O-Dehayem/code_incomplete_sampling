#' This function modifies a dataset of simulated islands by removing some entire random lineages.

 
#' Load necessary packages
library(ape)
library(DAISIE)
library(DDD)



remove_entire_lineages_randomly <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G) {
  data_sims_G19 <- list()  # Initialize an empty list to store the modified islands
  simulated_dataset_to_change_data <- data_sims_G  # Dataset to be modified
  simulated_tree_to_change_data <- data_sims_trees  # Trees to be modified

  one_island <- simulated_dataset_to_change_data[[i]]  # Extract data for one island
  one_island_tree <- simulated_tree_to_change_data[[i]]  # Extract trees for one island
  
  # Calculate the number of species to remove
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
  
  one_island1 <- one_island  # Create a copy of the island data to modify
  species_already_removed <- 0  # Initialize counter for removed species

  # Loop until the specified number of species has been removed
  while (species_already_removed < spec_to_remove) {
    for (j in 2:length(one_island)) {
      # Remove species until the required number is reached or only two branching times are left
      while (species_already_removed < spec_to_remove && length(one_island1[[j]]$branching_times) > 2) {
        trees <- one_island_tree[[j]]
        branching_times <- one_island1[[j]]$branching_times[-c(1, 2)]
        if (length(branching_times) > 1) {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          drop_result <- drop_youngest_species(trees, branching_times)
          brts <- drop_result$sorted_brts
          trees <- drop_result$tree_reduced
          New_branching_times <- c(one_island1[[j]]$branching_times[c(1, 2)], brts)
          one_island1[[j]]$branching_times <- New_branching_times
          branching_times <- c(one_island1[[j]]$branching_times[c(1, 2)], brts)
          species_already_removed <- species_already_removed + ((length_old_branching_times) - length(New_branching_times))
          one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + ((length_old_branching_times) - length(New_branching_times))
          one_island_tree[[j]] <- trees
        } else {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          New_branching_times <- one_island1[[j]]$branching_times[c(1, 2)]
          one_island1[[j]]$branching_times <- New_branching_times
          branching_times <- one_island1[[j]]$branching_times[c(1, 2)]
          species_already_removed <- species_already_removed + 1
          one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + 1
          one_island_tree[[j]] <- trees
        }
      }
      
      # Handle cases with exactly two branching times
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max) {
        if (one_island[[j]]$stac == 2) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        } else if (one_island[[j]]$stac == 3) {
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          one_island1[[j]]$stac <- 7
        }
      } else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max) {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max) {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
        species_already_removed <- species_already_removed + 2
        if (one_island[[j]]$stac == 2) {
          one_island1[[j]]$stac <- 5
        } else if (one_island[[j]]$stac == 3) {
          one_island1[[j]]$stac <- 7
        } else if (one_island[[j]]$stac == 4) {
          one_island1[[j]]$stac <- 1
        }
      }

      data_sims_G19 <- one_island1
    }
  }
  return(data_sims_G19)
}


