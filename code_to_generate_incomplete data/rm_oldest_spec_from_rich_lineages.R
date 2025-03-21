#' This function modifies a list of simulated datasets by removing the oldest species
#' from lineages that exceed a specified richness threshold.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages.
#' @param percentage_spec_to_remove Percentage of species to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
 


rm_oldest_spec_from_rich_linages <- function(all_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G) {
  data_sims_G1 <- list()  # Initialize an empty list to store the modified islands
  simulated_dataset_to_change_data <- data_sims_G  # Dataset to be modified
  simulated_tree_to_change_data <- all_trees  # Trees to be modified

  one_island <- simulated_dataset_to_change_data[[i]]  # Extract data for one island
  one_island_tree <- simulated_tree_to_change_data[[i]]  # Extract trees for one island
  
  # Calculate the number of species to remove
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
  
  one_island1 <- one_island  # Create a copy of the island data to modify
  species_already_removed <- 0  # Initialize counter for removed species

  # Loop until the specified number of species has been removed
  while (species_already_removed < spec_to_remove) {
    for (j in 2:length(one_island1)) {
      if (length(one_island[[j]]$branching_times) > rich_lineage_treshold && species_already_removed < spec_to_remove) {
        if (length(one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove) {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if (length(branching_times) > 1) {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S, function(x, y) which(y == x), y = trees$edge[, 2])], trees$tip.label)
            
            branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
            if (length(branch_to_drop_max) > 1) {
              branch_to_drop_max <- branch_to_drop_max[1]
            }
            tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
            tip_lengths <- tip_lengths[-(which(tip_lengths == max(tip_lengths)))[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        } else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max) {
          if (one_island[[j]]$stac == 2) {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          } else if (one_island[[j]]$stac == 3) {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
        } else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
    }
    data_sims_G1 <- one_island1
  }

  return(data_sims_G1)
}


}

rm_oldest_spec_from_rich_linages <- (data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
