# Two helper functions to impute taxonomic ranks into the taxonomic hierarchies returned by ncbi.


# function to return the taxonomic ranks straddling a particular 'main' taxonomic rank
get_straddling_tax_levels <- function(tax_level) {
  if(tax_level == 'genus') {
    return(c('subgenus', 'subtribe', 'tribe', 'subfamily'))
  } else if(tax_level == 'family') {
    return(c('subtribe', 'tribe', 'subfamily', 'superfamily', 'parvorder', 'infraorder', 'suborder'))
  } else if(tax_level == 'order') {
    return(c('superfamily', 'parvorder', 'infraorder', 'suborder', 'superorder', 'infraclass', 'subclass'))
  } else if(tax_level == 'class') {
    return(c('superorder', 'infraclass', 'subclass', 'superclass', 'subphylum'))
  } else if(tax_level == 'phylum') {
    return(c('superclass', 'subphylum', 'subkingdom'))
  }
}


# function to replace missing taxonomic rank: takes in species and missing taxonomic rank, return replacement name
fill_in_tax_level <- function(data, species, tax_level_missing){
  
  row.in.data <- which(data$sp.query == species) #get taxonomic data for species from main tax table
  
  # check the taxonomic info available for possible replacements
  tax_levels_to_use_for_replacement <- get_straddling_tax_levels(tax_level = tax_level_missing)
  possible_replacements <- as.character(data[row.in.data, tax_levels_to_use_for_replacement])
  possible_replacements <- possible_replacements[!is.na(possible_replacements)]
  
  if(length(possible_replacements) == 0) { # if there no available replacement, return NA
    replacement <- NA
  } else {
    
    # get the counts of taxonomic names across all host data
    counts <- as.numeric(common.tax.groups[
      match(possible_replacements, common.tax.groups$tax.name),]$count)
    
    # replacement is the name with biggest count
    replacement <- possible_replacements[which(counts == max(counts))]
    
    
    # if there is a tie for replacement 'commonness', then simply take the first one
    if(length(replacement) > 1) {
      replacement <- replacement[1]
    }
  }
  
  return(replacement)
}
