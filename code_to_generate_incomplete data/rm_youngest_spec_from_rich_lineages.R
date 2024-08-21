#' This function modifies a list of simulated datasets by removing the yungest species
#' from lineages that exceed a specified richness threshold.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages.
#' @param percentage_spec_to_remove Percentage of species to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
 
#' Load necessary packages
library(ape)
library(DAISIE)
library(DDD)



remove_youngest_spec_from_rich_lineages <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  # Initialize the output list
  data_sims_G2 <- list()
  
  # Reference the input datasets and trees
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees

  # Assuming that we need to iterate over each island dataset (index i should be defined elsewhere)
  for (i in seq_along(simulated_dataset_to_change_data)) {
    one_island <- simulated_dataset_to_change_data[[i]]
    one_island_tree <- simulated_tree_to_change_data[[i]]
    
    # Calculate the number of species to remove based on the given percentage
    spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
    
    # Create a copy of the dataset to modify
    one_island1 <- one_island
    species_already_removed <- 0
    
    # Continue modifying the dataset until the required number of species are removed
    while (species_already_removed < spec_to_remove) {
      for (j in 2:length(one_island1)) {
        
        # Case 1: Lineages with branching times greater than the richness threshold
        if (length(one_island[[j]]$branching_times) > rich_lineage_treshold && species_already_removed < spec_to_remove) {
          
          # Case 1.1: Lineages with more than two branching times
          if (length(one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove) {
            trees <- one_island_tree[[j]]
            branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
            
            # Drop a species with the minimum tip length from the tree
            if (length(branching_times) > 1) {
              length_old_branching_times <- length(one_island1[[j]]$branching_times)
              S <- length(branching_times) + 1
              tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                               function(x, y) which(y == x), y = trees$edge[,2])], trees$tip.label)
              
              branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
              if (length(branch_to_drop_min) > 1) {
                branch_to_drop_min <- branch_to_drop_min[1]
              }
              tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
              tip_lengths <- tip_lengths[-(which(tip_lengths == min(tip_lengths)))[1]]
              trees <- tree_reduced
              bt <- branching.times(tree_reduced)
              brts <- sort(as.vector(bt), decreasing = TRUE)
              New_branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
              one_island1[[j]]$branching_times <- New_branching_times
              species_already_removed <- species_already_removed + ((length_old_branching_times) - length(New_branching_times))
              one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + ((length_old_branching_times) - length(New_branching_times))
              one_island_tree[[j]] <- trees
            }
            else {
              # If only two branching times are present, remove one species
              length_old_branching_times <- length(one_island1[[j]]$branching_times)
              New_branching_times <- one_island1[[j]]$branching_times[c(1,2)]
              one_island1[[j]]$branching_times <- New_branching_times
              species_already_removed <- species_already_removed + 1
              one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + 1
              one_island_tree[[j]] <- trees
            }
          }
          
          # Case 1.2: Lineages with exactly two branching times, rich lineage, and not yet maxed out
          else if (length(one_island1[[j]]$branching_times) == 2 &&
                   length(one_island[[j]]$branching_times) > 2 &&
                   species_already_removed < spec_to_remove &&
                   one_island1[[j]]$branching_times[2] != col_max)
          {
            # Update stac and branching times based on original stac
            if(one_island[[j]]$stac == 2) {
              one_island1[[j]]$branching_times <- c(island_age, col_max)
              species_already_removed <- species_already_removed + 2
              one_island1[[j]]$stac <- 6
            } else if(one_island[[j]]$stac == 3) {
              one_island1[[j]]$branching_times <- c(island_age, col_max)
              species_already_removed <- species_already_removed + 2
              one_island1[[j]]$stac <- 7
            }
          }
          
          # Case 1.3: Lineages with exactly two branching times and maximum column value
          else if (length(one_island1[[j]]$branching_times) == 2 &&
                   length(one_island[[j]]$branching_times) > 2 &&
                   species_already_removed < spec_to_remove &&
                   one_island1[[j]]$branching_times[2] == col_max)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
          }
        }
      }
      # Update the list with modified dataset for the current island
      data_sims_G2[[i]] <- one_island1
    }
  }
  
  return(data_sims_G2)
}

# Example usage:
remove_youngest_spec_from_rich_lineages(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)
