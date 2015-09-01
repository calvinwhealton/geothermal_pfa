# function to assign drilling fluids into condensed groups

assignCondFluid <- function(df      #data frame with drilling fluid listed, fluid is the variable of the fluid name
                            ,flList # list of fluids to condense (group name assigned to each fluid)
                            ){
  # initializng variable to hold the condensed drilling fluid
  df$condFluid <- NA
  
  # looping over the groups in the data frame
  for(i in 1:length(names(flList))){
    
    # finding fluid names assigned to each group
    names_group <- flList[[i]]
    
    # looping over the assigned fluid names
    for(j in 1:as.numeric(summary(flList)[[i]])){
      
      # assigning group name when the fluid name matches a value in the list
      df$condFluid[which(df$fluid %in% names_group[j])] <- names(flList)[i]
    }
  }
  
  return(df)
}