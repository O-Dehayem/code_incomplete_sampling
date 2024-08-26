#' This function modifies a list of simulated datasets by removing random species
#' from lineages that fall below a specified richness threshold.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages.
#' @param percentage_spec_to_remove Percentage of species to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
#'
 
#' Load necessary packages
library(ape)
library(DAISIE)
library(DDD)




# Define the function remove_random_spec_from_poor_lineages
remove_random_spec_from_poor_lineages <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  # Initialize the output list
  data_sims_G6 <- list()
  
  # Assign the input data to local variables for clarity
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees

  # Select the first island dataset and tree (assuming i is predefined)
  one_island <- simulated_dataset_to_change_data[[i]]
  one_island_tree <- simulated_tree_to_change_data[[i]]
  
  # Calculate the number of species to remove based on the specified percentage
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum) * percentage_spec_to_remove, 0)
  
  # Create a copy of the original island dataset for modification
  one_island1 <- one_island
  species_already_removed <- 0
  
  # Continue the process until the desired number of species are removed
  while (species_already_removed < spec_to_remove)
  {
    # Iterate over each lineage in the island dataset starting from the second element
    for (j in 2:length(one_island1))
    {
      # Check if the lineage has fewer branching times than the threshold and species need to be removed
      if (length(one_island[[j]]$branching_times) <= rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        # If the branching times are exactly 2 and meet specific conditions, remove species
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && 
            species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          # Update branching times and status
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          
          # Update the stac value based on the original dataset
          if (one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$stac <- 5
          }
          else if (one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$stac <- 7
          }
          else if (one_island[[j]]$stac == 4)
          {
            one_island1[[j]]$stac <- 1
          }
        }
        # If branching times are greater or equal to 3, randomly remove a species
        else if (length(one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          
          if (length(branching_times) > 1)
          {
            # Calculate the number of tips and their lengths
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S, function(x,y) which(y == x), y = trees$edge[,2])], trees$tip.label)

            # Randomly select a tip length to drop
            set.seed(5)
            sampled_length <- sample(tip_lengths, 1)

            # Get the indices of elements equal to the sampled length
            matching_indices <- which(tip_lengths == sampled_length)

            # Drop the species
            tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices])[1])
            tip_lengths <- tip_lengths[-matching_indices[1]]
            trees <- tree_reduced
            
            # Recalculate the branching times
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
            New_branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
            one_island1[[j]]$branching_times <- New_branching_times
            branching_times <- c(one_island1[[j]]$branching_times[c(1,2)], brts)
            
            # Update species removal count and status
            species_already_removed <- species_already_removed + (length_old_branching_times - length(New_branching_times))
            one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + (length_old_branching_times - length(New_branching_times))
            one_island_tree[[j]] <- trees
          }
          else
          {
            # Handle cases with fewer branching times
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            New_branching_times <- one_island1[[j]]$branching_times[c(1,2)]
            one_island1[[j]]$branching_times <- New_branching_times
            branching_times <- one_island1[[j]]$branching_times[c(1,2)]
            species_already_removed <- species_already_removed + 1
            one_island1[[j]]$missing_species <- one_island1[[j]]$missing_species + 1
            one_island_tree[[j]] <- trees
          }
        }
        # Handle cases where branching times are 2, but the original has more branching times
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && 
                 species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          if (one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if (one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
        }
        # Handle cases where branching times are 2 and meet specific conditions
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && 
                 one_island1[[j]]$branching_times[2] == col_max)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
    }
    # Update the output dataset
    data_sims_G6 <- one_island1
  }
  
  # Return the modified dataset
  return (data_sims_G6)
}

# Execute the function with the provided parameters
remove_random_spec_from_poor_lineages(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)
