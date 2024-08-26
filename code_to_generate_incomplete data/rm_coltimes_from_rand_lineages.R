#' Remove Colonization Times from Random Lineages in Simulated Datasets
#'
#' This function modifies a list of simulated datasets by removing outgroup species
#' from random lineages based on specified conditions.
#'
#' @param rich_lineage_threshold The threshold for richness of lineages (not directly used in this version).
#' @param percentage_coltimes_to_remove Percentage of colonization times to remove based on total species.
#' @param data_sims_G A list of simulated datasets.
#' @return A modified list of simulated datasets with adjusted colonization times.
#'

#' Load necessary library
library(ape)
library(DAISIE)
library(DDD)



remove_coltimes_from_random <- function(rich_lineage_threshold, percentage_coltimes_to_remove, data_sims_G) {

  # Initialize an empty list to store modified datasets
  data_sims_G18 <- list()

  # Iterate over each dataset in data_sims_G
  for (i in seq_along(data_sims_G)) {
    # Make a copy of the current dataset
    one_island <- data_sims_G[[i]]

    # Calculate number of colonization times to remove
    total_species <- sum(Number_of_species(data_sims_G)[[i]])
    coltimes_to_remove <- round(total_species * percentage_coltimes_to_remove, 0)

    # Initialize counter for number of colonization times already removed
    number_of_coltimes_already_removed <- 0

    # Loop until the desired number of colonization times are removed
    while (number_of_coltimes_already_removed < coltimes_to_remove) {
      # Iterate over each lineage in the current dataset
      for (j in 2:length(one_island)) {
        # Check if enough colonization times have been removed
        if (number_of_coltimes_already_removed < coltimes_to_remove) {

          # Remove outgroup species based on lineage status and branching times
          if (length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 4 && one_island[[j]]$branching_times[2] != col_max) {
            # Adjust colonization time and update lineage status
            one_island[[j]]$branching_times[2] <- col_max
            number_of_coltimes_already_removed <- number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 1
          } else if (one_island[[j]]$stac == 3 && one_island[[j]]$branching_times[2] != col_max) {
            # Adjust colonization time and update lineage status
            one_island[[j]]$branching_times[2] <- col_max
            number_of_coltimes_already_removed <- number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 7
          } else if (length(one_island[[j]]$branching_times) > 2 && one_island[[j]]$stac == 2 && one_island[[j]]$branching_times[2] != col_max) {
            # Adjust colonization time and update lineage status
            one_island[[j]]$branching_times[2] <- col_max
            number_of_coltimes_already_removed <- number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 6
          } else if (length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 2 ) {
            # Adjust colonization time and update lineage status
            one_island[[j]]$branching_times[2] <- col_max
            number_of_coltimes_already_removed <- number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 5
          }

        }
      }
    }

    # Store the modified dataset in the list of modified datasets
    data_sims_G18[[i]] <- one_island
  }

  # Return the list of modified datasets
  return(data_sims_G18)
}

remove_coltimes_from_random <- data_sims_G18(rich_lineage_treshold,coltimes_to_remove,data_sims_G)

