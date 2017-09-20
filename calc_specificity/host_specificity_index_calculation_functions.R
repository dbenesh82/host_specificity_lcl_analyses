
#functions for calculating host specificity index

td.level <- function(taxa){
  #takes in a vector of taxonomic information, e.g. a list of genera
  #creates all pairwise combinations in tax list and differences as logical vector
  #helper function for hs.index function below
  x <- combn(taxa, 2)
  td <- x[1,] != x[2,]
  return(td)
}


hs.index <- function(species, genus, family, order, class, phylum) {
  #takes in taxonomic information for a list of species
  #returns the index of host specificity proposed by Poulin and Mouillot 2003
  #also returns the var of index
  if(length(species)==1) { #if just a single species, return value of 1
    hs <- 1
    var.hs <- NA
  } else { #calculate index
    num.sp <- length(species) #number of species
    td <- td.level(species) + td.level (genus) + td.level(family) + td.level(order) +
      td.level(class) + td.level(phylum)
    sum.td <- sum(td)
    hs <- (2*sum.td/(num.sp * (num.sp-1))) #index
    var.hs <- 2*sum((td - mean(td))^2) / (num.sp * (num.sp-1)) #variance of index; gives right result, but had to double it
  }
  return(c(hs, var.hs))
}


test.case<-data.frame(species = letters[1:4], 
                      genus = c("a","a","b","b"),
                      family = c("a","a","a","a"),
                      order = rep("a", 4),
                      class = rep("a", 4),
                      phylum = rep("a", 4))

with(test.case, hs.index(species, genus, family, order, class, phylum))


test.case<-data.frame(species = letters[1:7], 
                      genus = c("a","a","a","b","b","b","c"),
                      family = c(rep("a", 6), "b"),
                      order = rep("a", 7),
                      class = rep("a", 7),
                      phylum = rep("a", 7))

with(test.case, hs.index(species, genus, family, order, class, phylum))

#looks like it works; returns same values as in Poulin and Mouillot paper
rm(test.case)

