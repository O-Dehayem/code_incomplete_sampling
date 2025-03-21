# Load necessary library



### function to convert branching times in datalist to trees with fixes topology
# Function to convert branching times in data_list to trees with fixed topology
convert_branching_times_to_trees <- function(data_list) {
  
  # Initialize an empty list to store the resulting trees
  data_list_trees <- list()
  
  # Loop through each element in the dataset starting from the second element
  for (j in 2:length(data_list)) {
    
    # Get the branching times for the current dataset
    branching_times <- data_list[[j]]$branching_times
    
    # Initialize a variable to store the current tree
    one_tree <- NULL
    
    # Check the length of branching times and convert accordingly
    if (length(branching_times) > 3) {
      # Convert branching times to a phylogenetic tree, excluding the first two branching times
      one_tree <- DDD::brts2phylo(branching_times[-c(1, 2)])
    } else if (length(branching_times) == 3) {
      # Convert branching times to a phylogenetic tree, excluding the first branching time and setting root = TRUE
      one_tree <- DDD::brts2phylo(branching_times[-1], root = TRUE)
    } else if (length(branching_times) == 2) {
      # Convert branching times to a phylogenetic tree, excluding the first branching time
      one_tree <- DDD::brts2phylo(branching_times[-1])
    }
    
    # Store the resulting tree in the list of trees
    data_list_trees[[j]] <- one_tree
  }
  
  # Return the list of trees
  return(data_list_trees)
}

# Example usage
data_list_trees <- convert_branching_times_to_trees(data_list)
data_list_trees
