#' This function modifies a list of simulated datasets by removing some entire lineages
#' from lineages that fall below a specified richness threshold.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages.
#' @param percentage_spec_to_remove Percentage of species to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
#'
 


remove_entire_rich_lineages <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  # Initialize an empty list to store the modified island datasets
  data_sims_G14 <- list()
  
  # Assign the provided datasets to local variables for manipulation
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees

  # Assuming 'i' is defined somewhere else in the scope
  one_island <- simulated_dataset_to_change_data[[i]]
  
  # Create 'one_island_tree', an island dataset with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  
  # Calculate the number of species to remove
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
  
  # Initialize a variable to store the modified dataset
  one_island1 <- one_island
  species_already_removed <- 0
  
  # Loop until the desired number of species are removed
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island))
    {
      if (length(one_island[[j]]$branching_times) > rich_lineage_treshold)
      {
        # Remove species while the number of species to remove is not met and branching times are greater than 2
        while (species_already_removed < spec_to_remove && length(one_island1[[j]]$branching_times) > 2)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          
          if (length(branching_times) > 1)
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            drop_result <- drop_youngest_species(trees, branching_times)
            brts <- drop_result$sorted_brts
            trees <- drop_result$tree_reduced
            New_branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
            one_island1[[j]]$branching_times <- New_branching_times
            branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
            species_already_removed <- species_already_removed + ((length_old_branching_times) - length(New_branching_times))
            one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + ((length_old_branching_times) - length(New_branching_times))
            one_island_tree[[j]] <- trees
          }
          else
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            New_branching_times <- one_island1[[j]]$branching_times[c(1,2)]
            one_island1[[j]]$branching_times <- New_branching_times
            branching_times <- one_island1[[j]]$branching_times[c(1,2)]
            species_already_removed <- species_already_removed + 1
            one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + 1
            one_island_tree[[j]] <- trees
          }
        }
        
        # Handle specific cases based on the branching times and stac value
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          if (one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if (one_island[[j]]$stac == 3)
          {
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            one_island1[[j]]$stac <- 7
          }
        }
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
      
      # Update the final modified dataset list
      data_sims_G14 <- one_island1
    }
  }
  return (data_sims_G14)
}
# Execute the function with the given parameters
remove_entire_rich_lineages(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

