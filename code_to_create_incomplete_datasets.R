#packages needed
library(ape)



data_sims_trees <- list()
for ( i in seq_along(data_sims_G))
{
  trees <- list()
  one_dataset <- data_sims_G[[i]]
  for (j in 2:length(one_dataset)) 
  {
    
    if (length(one_dataset[[j]]$branching_times) > 3)
    {
      one_tree <- DDD::brts2phylo(one_dataset[[j]]$branching_times[-c(1,2)])   
    }
    else if  (length(one_dataset[[j]]$branching_times) == 3) 
    {
      one_tree <- DDD::brts2phylo(one_dataset[[j]]$branching_times[-1], root = TRUE)   
    }
    else if  (length(one_dataset[[j]]$branching_times) == 2) 
    {
      one_tree <- DDD::brts2phylo(one_dataset[[j]]$branching_times[-1])   
    }
    
    trees[[j]] <- one_tree
  }
  data_sims_trees[[i]] <- trees
}

###############
drop_youngest_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
  # Check if branch_to_drop_min has multiple values
  if (length(branch_to_drop_min) > 1) {
    #Drop only one of the youngest species
    branch_to_drop_min <- branch_to_drop_min[1]
  }
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}

###############
drop_oldest_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
  # Check if branch_to_drop_min has multiple values
  if (length(branch_to_drop_max) > 1) {
    # Drop only one of the youngest species
    branch_to_drop_max <- branch_to_drop_max[1]
  }
  tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}

###############
drop_random_species <- function(trees, branching_times)
{
  
  S <- length(branching_times) + 1
  tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                   function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
  
  sampled_length <- sample(tip_lengths, 1)
  
  # Get the indices of elements equal to the sampled length
  matching_indices <- which(tip_lengths == sampled_length)
  
  # drop the species
  
  
  tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices]))
  brts <- branching.times(tree_reduced)
  sorted_brts <- sort(as.vector(brts), decreasing = TRUE)
  return(list(tip_lengths = tip_lengths, tree_reduced = tree_reduced, sorted_brts = sorted_brts))
  
}


rich_lineage_treshold <- 4 
percentage_spec_to_remove <- 0.4
island_age <- 4
col_max <- island_age - 1e-05


Branching_times <-  function(dataset)
{
  Branching_times_list = list()
  
  for (i in 1:length(dataset))
  {
    one_island <- dataset[[i]]
    Branching_times_one_island <- data.frame()
    for (j in 2:length(one_island))
    {
      data_Bt = data.frame(length(one_island[[j]]$branching_times))
      Branching_times_one_island = rbind(Branching_times_one_island, data_Bt)
    }
    Branching_times_list[i] <- Branching_times_one_island
    Number_of_col_vector <- sapply(dataset,length) - 1
  }
  Branching_times_list
}


Number_of_species <-  function(dataset)
{
  Number_of_species_list = list()
  
  for (i in 1:length(dataset))
  {
    one_island <- dataset[[i]]
    Branching_times_one_island <- data.frame()
    for (j in 2:length(one_island))
    {
      data_Bt = data.frame(length(one_island[[j]]$branching_times) - 1)
      Branching_times_one_island = rbind(Branching_times_one_island, data_Bt)
    }
    Number_of_species_list[i] <- Branching_times_one_island
  }
  Number_of_species_list
}




rich_lineage_treshold <- 4 
percentage_spec_to_remove <- 0.6






data_sims_G1 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G1 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      if (length (one_island[[j]]$branching_times) > rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        
        if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if (length(branching_times) > 1)
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            
            branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
            if (length(branch_to_drop_max) > 1) 
            {
              branch_to_drop_max <- branch_to_drop_max[1]
            }
            tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
            tip_lengths <- tip_lengths [-(which(tip_lengths == max(tip_lengths)))[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2  && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times = c(island_age, col_max)
            species_already_removed = species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
        else if (length(one_island1[[j]]$branching_times) == 2  && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max)
        {
          
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
    }
    data_sims_G1 <- one_island1
  }
  
  return (data_sims_G1)
  
}

data_sims_G113 <- data_sims_G1(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G2 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G2 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island 
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      if (length (one_island[[j]]$branching_times) > rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        
        if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if (length( branching_times) > 1 )
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            
            branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
            if (length(branch_to_drop_min) > 1) 
            {
              branch_to_drop_min <- branch_to_drop_min[1]
            }
            tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
            tip_lengths <- tip_lengths [-(which(tip_lengths == min(tip_lengths)))[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2  && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times = c(island_age, col_max)
            species_already_removed = species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
        else if (length(one_island1[[j]]$branching_times) == 2  && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
    }
    data_sims_G2 <- one_island1
  }
  
  return (data_sims_G2)
  
}

data_sims_G222 <- data_sims_G2(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G3 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G3 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      if (length (one_island[[j]]$branching_times) > rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        
        if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if (length( branching_times) > 1 )
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            set.seed (5)
            sampled_length <- sample(tip_lengths, 1)
            
            # Get the indices of elements equal to the sampled length
            matching_indices <- which(tip_lengths == sampled_length)
            
            # drop the species
            
            
            tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices])[1])
            tip_lengths <- tip_lengths [-matching_indices[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
    }
    data_sims_G3 <- one_island1
  }
  
  return (data_sims_G3)
  
}

data_sims_G33 <- data_sims_G3(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G4 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G4 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length (one_island[[j]]$branching_times) <= rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          one_island1[[j]]$branching_times = c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          if(one_island[[j]]$stac == 2) 
          {
            one_island1[[j]]$stac <- 5
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$stac <- 7
          }
          else if(one_island[[j]]$stac == 4)
          {
            one_island1[[j]]$stac <- 1
          }
          
        }
        
        else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if( length( branching_times) > 1 )
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            
            branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
            if (length(branch_to_drop_max) > 1) 
            {
              branch_to_drop_max <- branch_to_drop_max[1]
            }
            tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
            tip_lengths <- tip_lengths [-(which(tip_lengths == max(tip_lengths)))[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3) {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
      }
      
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
      
      
    }
    data_sims_G4 <- one_island1
  }
  
  return (data_sims_G4)
}

data_sims_G4 <- data_sims_G4(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G5 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G5 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length (one_island[[j]]$branching_times) <= rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          if(one_island[[j]]$stac == 2) 
          {
            one_island1[[j]]$stac <- 5
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$stac <- 7
          }
          else if(one_island[[j]]$stac == 4)
          {
            one_island1[[j]]$stac <- 1
          }
          
        }
        
        else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if( length( branching_times) > 1 )
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            
            branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
            if (length(branch_to_drop_min) > 1) 
            {
              branch_to_drop_min <- branch_to_drop_min[1]
            }
            tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
            tip_lengths <- tip_lengths [-(which(tip_lengths == min(tip_lengths)))[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3) {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
      
      
      
    }
    data_sims_G5 <- one_island1
  }
  
  return (data_sims_G5)
}


data_sims_G5 <- data_sims_G5(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G6 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G6 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length (one_island[[j]]$branching_times) <= rich_lineage_treshold && species_already_removed < spec_to_remove)
      {
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          one_island1[[j]]$branching_times = c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          if(one_island[[j]]$stac == 2) 
          {
            one_island1[[j]]$stac <- 5
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$stac <- 7
          }
          else if(one_island[[j]]$stac == 4)
          {
            one_island1[[j]]$stac <- 1
          }
          
        }
        
        else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
        {
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if( length( branching_times) > 1 )
          {
            length_old_branching_times <- length(one_island1[[j]]$branching_times)
            S <- length(branching_times) + 1
            tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                             function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
            
            set.seed (5)
            sampled_length <- sample(tip_lengths, 1)
            
            # Get the indices of elements equal to the sampled length
            matching_indices <- which(tip_lengths == sampled_length)
            
            # drop the species
            
            
            tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices])[1])
            tip_lengths <- tip_lengths [-matching_indices[1]]
            trees <- tree_reduced
            bt <- branching.times(tree_reduced)
            brts <- sort(as.vector(bt), decreasing = TRUE)
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
        
        else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3) {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 7
          }
          
        }
        
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
        }
      }
      
      
      
    }
    data_sims_G6 <- one_island1
  }
  
  return (data_sims_G6)
}

data_sims_G6 <- data_sims_G6(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)




data_sims_G7 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G7 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        one_island1[[j]]$branching_times = c(island_age, col_max)
        species_already_removed <- species_already_removed + 2
        if(one_island[[j]]$stac == 2) 
        {
          one_island1[[j]]$stac <- 5
        }
        else if(one_island[[j]]$stac == 3)
        {
          one_island1[[j]]$stac <- 7
        }
        else if(one_island[[j]]$stac == 4)
        {
          one_island1[[j]]$stac <- 1
        }
      }
      
      else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
      {
        trees <- one_island_tree[[j]]
        branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
        if( length( branching_times) > 1 )
        {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          S <- length(branching_times) + 1
          tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                           function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
          
          
          branch_to_drop_max <- names(which(tip_lengths == max(tip_lengths)))
          if (length(branch_to_drop_max) > 1) 
          {
            branch_to_drop_max <- branch_to_drop_max[1]
          }
          tree_reduced <- ape::drop.tip(trees, branch_to_drop_max)
          tip_lengths <- tip_lengths [-(which(tip_lengths == max(tip_lengths)))[1]]
          trees <- tree_reduced
          bt <- branching.times(tree_reduced)
          brts <- sort(as.vector(bt), decreasing = TRUE)
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
      
      else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        if(one_island[[j]]$stac == 2)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        }
        else if(one_island[[j]]$stac == 3) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 7
        }
        
      }
      
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
      
      
      
      
    }
    data_sims_G7 <- one_island1
  }
  
  return (data_sims_G7)
}


data_sims_G7 <- data_sims_G7(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

data_sims_G8 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G8 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length(one_island1[[j]]$branching_times)==2 && length(one_island[[j]]$branching_times)==2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        one_island1[[j]]$branching_times = c(island_age, col_max)
        species_already_removed = species_already_removed + 2
        if(one_island[[j]]$stac == 2) 
        {
          one_island1[[j]]$stac <- 5
        }
        else if(one_island[[j]]$stac == 3)
        {
          one_island1[[j]]$stac <- 7
        }
        else if(one_island[[j]]$stac == 4)
        {
          one_island1[[j]]$stac <- 1
        }
        
      }
      
      else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
      {
        trees <- one_island_tree[[j]]
        branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
        if( length( branching_times) > 1 )
        {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          S <- length(branching_times) + 1
          tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                           function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
          
          
          branch_to_drop_min <- names(which(tip_lengths == min(tip_lengths)))
          if (length(branch_to_drop_min) > 1) 
          {
            branch_to_drop_min <- branch_to_drop_min[1]
          }
          tree_reduced <- ape::drop.tip(trees, branch_to_drop_min)
          tip_lengths <- tip_lengths [-(which(tip_lengths == min(tip_lengths)))[1]]
          trees <- tree_reduced
          bt <- branching.times(tree_reduced)
          brts <- sort(as.vector(bt), decreasing = TRUE)
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
      
      else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) >2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        if(one_island[[j]]$stac == 2)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        }
        else if(one_island[[j]]$stac == 3) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 7
        }
        
      }
      
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
      
      
    }
    data_sims_G8 <- one_island1
  }
  
  return (data_sims_G8)
}


data_sims_G8 <- data_sims_G8(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)


data_sims_G9 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  
  data_sims_G9 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island1))
    {
      
      if (length(one_island1[[j]]$branching_times)==2 && length(one_island[[j]]$branching_times)==2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
        species_already_removed = species_already_removed + 2
        if(one_island[[j]]$stac == 2) 
        {
          one_island1[[j]]$stac <- 5
        }
        else if(one_island[[j]]$stac == 3)
        {
          one_island1[[j]]$stac <- 7
        }
        else if(one_island[[j]]$stac == 4)
        {
          one_island1[[j]]$stac <- 1
        }
        
      }
      
      else if (length (one_island1[[j]]$branching_times) >= 3 && species_already_removed < spec_to_remove)
      {
        trees <- one_island_tree[[j]]
        branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
        if( length( branching_times) > 1 )
        {
          length_old_branching_times <- length(one_island1[[j]]$branching_times)
          S <- length(branching_times) + 1
          tip_lengths <- setNames(trees$edge.length[sapply(1:S,
                                                           function(x,y) which(y == x),y = trees$edge[,2])],trees$tip.label)
          
          set.seed (5)
          sampled_length <- sample(tip_lengths, 1)
          
          # Get the indices of elements equal to the sampled length
          matching_indices <- which(tip_lengths == sampled_length)
          
          # drop the species
          
          
          tree_reduced <- ape::drop.tip(trees, names(tip_lengths[matching_indices])[1])
          tip_lengths <- tip_lengths [-matching_indices[1]]
          trees <- tree_reduced
          bt <- branching.times(tree_reduced)
          brts <- sort(as.vector(bt), decreasing = TRUE)
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
      
      else if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        if(one_island[[j]]$stac == 2)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        }
        else if(one_island[[j]]$stac == 3) {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 7
        }
        
      }
      
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
      }
      
      
      
      
    }
    data_sims_G9 <- one_island1
  }
  
  return (data_sims_G9)
}


data_sims_G9 <- data_sims_G9(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)


data_sims_G14 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  data_sims_G14 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that will be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island))
    { 
      if (length (one_island[[j]]$branching_times) > rich_lineage_treshold)
      {
        #length_old_branching_times = length(one_island[[j]]$branching_times)
        
        while (species_already_removed < spec_to_remove && length (one_island1[[j]]$branching_times) > 2)
        {
          
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if( length( branching_times) > 1 )
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
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3) 
          {
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            one_island1[[j]]$stac <- 7
          }
        }
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          
        }
        
        
      }
      
      data_sims_G14 <- one_island1
    }
  }
  return (data_sims_G14)
}
data_sims_G14 <- data_sims_G14(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)



data_sims_G15 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  data_sims_G15 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- all_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that wil be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island))
    { 
      if (length (one_island[[j]]$branching_times) <= rich_lineage_treshold)
      {
        #length_old_branching_times = length(one_island[[j]]$branching_times)
        
        while (species_already_removed < spec_to_remove && length (one_island1[[j]]$branching_times) > 2)
        {
          
          trees <- one_island_tree[[j]]
          branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
          if( length( branching_times) > 1 )
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
        if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          
          if(one_island[[j]]$stac == 2)
          {
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3) 
          {
            species_already_removed <- species_already_removed + 2
            one_island1[[j]]$branching_times <- c(island_age, col_max)
            one_island1[[j]]$stac <- 7
          }
        }
        else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          
        }
        if ( length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times)==2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          if(one_island[[j]]$stac == 2) 
          {
            one_island1[[j]]$stac <- 5
          }
          else if(one_island[[j]]$stac == 3)
          {
            one_island1[[j]]$stac <- 7
          }
          else if(one_island[[j]]$stac == 4)
          {
            one_island1[[j]]$stac <- 1
          }
          
        }
      }
      
      data_sims_G15 <- one_island1
    }
  }
  return (data_sims_G15)
}
data_sims_G15 <- data_sims_G15(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)
#######

data_sims_G16 <-  function(rich_lineage_treshold, coltimes_to_remove, data_sims_G)
{
  data_sims_G16 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  
  for (i in 1:length(data_sims_G))
  {
    one_island <- simulated_dataset_to_change_data[[i]]
    coltimes_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_coltimes_to_remove, 0)
    number_of_coltimes_already_removed = 0
    while (number_of_coltimes_already_removed < coltimes_to_remove)
    {
      for (j in 2:length(one_island))
      {
        if (length (one_island[[j]]$branching_times) > rich_lineage_treshold && number_of_coltimes_already_removed < coltimes_to_remove )
        {
          
          if(one_island[[j]]$stac == 2 && one_island[[j]]$branching_times[2] != col_max)
          {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3 && one_island[[j]]$branching_times[2] != col_max) {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 7
          }
          
        }
        
      }
      
    }
    data_sims_G16[[i]] <- one_island
  }
  return (data_sims_G16)
}
data_sims_G16 <- data_sims_G16(rich_lineage_treshold,coltimes_to_remove, data_sims_G)


data_sims_G17 <- function(rich_lineage_treshold, coltimes_to_remove, data_sims_G)
{
  data_sims_G17 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  
  for (i in 1:length(data_sims_G))
  {
    one_island <- simulated_dataset_to_change_data[[i]]
    number_of_coltimes_already_removed = 0
    coltimes_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_coltimes_to_remove, 0)
    while (number_of_coltimes_already_removed < coltimes_to_remove)
    {
      for (j in 2:length(one_island))
      {
        if (length (one_island[[j]]$branching_times) <= rich_lineage_treshold && number_of_coltimes_already_removed < coltimes_to_remove)
        {
          
          if(length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 4 && one_island[[j]]$branching_times[2] != col_max)
          {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 1
          }
          else if(length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 2 && one_island[[j]]$branching_times[2] != col_max) {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 5
          }
          else if(length(one_island[[j]]$branching_times) > 2 && one_island[[j]]$stac == 2 && one_island[[j]]$branching_times[2] != col_max) {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 6
          }
          else if(one_island[[j]]$stac == 3 && one_island[[j]]$branching_times[2] != col_max) {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 7
          }
          
        }
        
      }
    }
    
    data_sims_G17[[i]] <- one_island
  }
  return (data_sims_G17)
}
data_sims_G17 <- data_sims_G17(rich_lineage_treshold,coltimes_to_remove,data_sims_G)


data_sims_G18 <- function(rich_lineage_treshold, coltimes_to_remove, data_sims_G)
{
  data_sims_G18 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  
  for (i in 1:length(data_sims_G))
  {
    one_island <- simulated_dataset_to_change_data[[i]]
    coltimes_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_coltimes_to_remove, 0)
    number_of_coltimes_already_removed = 0
    while (number_of_coltimes_already_removed < coltimes_to_remove)
    {
      for (j in 2:length(one_island))
      {
        if (number_of_coltimes_already_removed < coltimes_to_remove)
        {
          
          if(length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 4 && one_island[[j]]$branching_times[2] != col_max)
          {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 1
          }
          else if(one_island[[j]]$stac == 3 && one_island[[j]]$branching_times[2] != col_max)
          {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 7
          }
          else if(length(one_island[[j]]$branching_times) > 2 && one_island[[j]]$stac == 2 && one_island[[j]]$branching_times[2] != col_max)
          {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 6
          }
          else if(length(one_island[[j]]$branching_times) == 2 && one_island[[j]]$stac == 2 ) {
            one_island[[j]]$branching_times[2] = col_max
            number_of_coltimes_already_removed = number_of_coltimes_already_removed + 1
            one_island[[j]]$stac <- 5
          }
          
        } 
      }
    }
    
    data_sims_G18[[i]] <- one_island
  }
  return (data_sims_G18)
}
data_sims_G18 <- data_sims_G18(rich_lineage_treshold,coltimes_to_remove,data_sims_G)




data_sims_G19 <- function(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, data_sims_G)
{
  data_sims_G19 <- list()
  simulated_dataset_to_change_data <- data_sims_G
  simulated_tree_to_change_data <- data_sims_trees
  
  one_island <- simulated_dataset_to_change_data[[i]]
  # create one_island_tree, an island data set with trees instead of branching times
  one_island_tree <- simulated_tree_to_change_data[[i]]
  spec_to_remove <- round(sapply(Number_of_species(data_sims_G)[i], sum)*percentage_spec_to_remove, 0)
  # create another island data set one_island1 that wil be modify
  one_island1 <- one_island
  species_already_removed <- 0
  while (species_already_removed < spec_to_remove)
  {
    for (j in 2:length(one_island))
    { 
      
      while (species_already_removed < spec_to_remove && length (one_island1[[j]]$branching_times) > 2)
      {
        
        trees <- one_island_tree[[j]]
        branching_times <- one_island1[[j]]$branching_times[-c(1,2)]
        if( length( branching_times) > 1 )
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
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) > 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        
        if(one_island[[j]]$stac == 2)
        {
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$stac <- 6
        }
        else if(one_island[[j]]$stac == 3) 
        {
          species_already_removed <- species_already_removed + 2
          one_island1[[j]]$branching_times <- c(island_age, col_max)
          one_island1[[j]]$stac <- 7
        }
      }
      else if (length(one_island1[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] == col_max )
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
        
      }
      if (length(one_island1[[j]]$branching_times) == 2 && length(one_island[[j]]$branching_times) == 2 && species_already_removed < spec_to_remove && one_island1[[j]]$branching_times[2] != col_max)
      {
        one_island1[[j]]$branching_times <- c(island_age, col_max)
        species_already_removed = species_already_removed + 2
        if(one_island[[j]]$stac == 2) 
        {
          one_island1[[j]]$stac <- 5
        }
        else if(one_island[[j]]$stac == 3)
        {
          one_island1[[j]]$stac <- 7
        }
        else if(one_island[[j]]$stac == 4)
        {
          one_island1[[j]]$stac <- 1
        }
        
      }
      
      
      data_sims_G19 <- one_island1
    }
  }
  return (data_sims_G19)
}
data_sims_G19 <- data_sims_G19(data_sims_trees, rich_lineage_treshold, percentage_spec_to_remove, Galapagos_sims_G)

##############
