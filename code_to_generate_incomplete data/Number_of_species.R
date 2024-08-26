#' The function takes a dataset as input and returns a list where each element corresponds to an island in the dataset.
#' For each island, the function calculates the number of species by counting the number of branching times minus one for 
#' lineages.



Number_of_species <- function(dataset) {
  # Initialize an empty list to store the number of species for each island
  Number_of_species_list <- list()

  # Loop over each island in the dataset
  for (i in 1:length(dataset)) {
    # Extract the data for one island
    one_island <- dataset[[i]]
    
    # Initialize a data frame to store the number of species for each group within the island
    Branching_times_one_island <- data.frame()
    
    # Loop over each group of species/trees on the island, starting from the second element
    for (j in 2:length(one_island)) {
      # Calculate the number of species for the current group by subtracting one from the number of branching times
      data_Bt <- data.frame(length(one_island[[j]]$branching_times) - 1)
      
      # Add the result to the data frame for the current island
      Branching_times_one_island <- rbind(Branching_times_one_island, data_Bt)
    }
    
    # Add the data frame for the current island to the list of results
    Number_of_species_list[[i]] <- Branching_times_one_island
  }
  
  # Return the list of results
  return(Number_of_species_list)
}
