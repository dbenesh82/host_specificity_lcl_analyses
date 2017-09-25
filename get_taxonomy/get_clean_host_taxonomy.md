Host specificity analyses - get and clean the host taxonomy
================

**Background**: I use a [database](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1680/suppinfo) of helminth (parasitic worm) life cycles to examine the diversity of hosts exploited by these parasites. Complex life cycle parasites infect dissimilar hosts over the course of their life cycles, e.g. first an invertebrate and then a vertebrate. This generalism is presumed to be costly. This script queries a taxonomic database to retrieve taxonomic information for the hosts in the database. This information is then used to calculate host specificity indexes.

First, import the libraries and the hosts in the life cycle database.

``` r
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(taxize)

options(stringsAsFactors = FALSE)

dataH <- read.csv(file="../data/CLC_database_hosts.csv", header = TRUE, sep=",")
```

Then we isolate the set of host species and genera, for which we desire taxonomies.

``` r
# extract genus from host species name
dataH <- mutate(dataH, h.genus = unlist(strsplit(as.character(Host.species), split=" .*$")))%>%
  select(Host.species, h.genus)

# get the set of host species and host genera
h.sp.vec <- unique(dataH$Host.species)
h.sp.vec <- sort(na.omit(h.sp.vec))

g.sp.vec <- unique(dataH$h.genus)
g.sp.vec <- sort(na.omit(g.sp.vec))
```

To get the taxonomies, we query the [NCBI database](https://www.ncbi.nlm.nih.gov/taxonomy) for each host species. This code is simply to get an idea how long it takes. It's blocked out because it does not need to be run everytime.

``` r
# QUERY NCBI DATABASE
# s1 <- system.time(ids.for.classif<-get_uid_(h.sp.vec[1:10]))
# s2 <- system.time(ids.for.classif<-get_uid_(h.sp.vec[1:50]))
# s3 <- system.time(ids.for.classif<-get_uid_(h.sp.vec[1:100]))
# s4 <- system.time(ids.for.classif<-get_uid_(h.sp.vec[1:500]))
# 
# tst <- data.frame(sp = c(10,50,100,500), time = c(s1['elapsed'], s2['elapsed'], s3['elapsed'], s4['elapsed']))
# ct <- coef(lm(time ~ sp, data=tst))
# 
# min <- (ct[1] + ct[2] * length(h.sp.vec) ) / 60
# 
# print(
#   paste("It's estimated to take", round(min, 1), "minutes to query all the species")
# )
# rm(s1, s2, s3, s4, tst, ct, min)
# on my laptop, it is predicted to take about an hour
```

For 3500+ host species, we query the NCBI database to get the uniques ids for each. It takes a while. The last time this script was run and the database was queried was 2017-09-25 15:45:00.

``` r
ids.for.classif <- get_uid_(h.sp.vec)
save.image("workspace_after_tax_downloads.RData")
```

This returns a list for each species that includes a 'uid', the unique identifier in NCBI. We want to extract uid, as well as some other data. The next block loops through each species to extract ids, the string queried (species or genus name), and the taxonomic level of the query (species or genus). When host species were not found in the NCBI database, the genus name was queried.

``` r
# empty character vectors for storing data collected in loop
ncbi.ids <- character() # unique ncbi id
tg.vec <- character() # the species name queried
q.rank <- character() # query at species or genus level

# loop through each species
for(i in seq_along(ids.for.classif)) {
  
  tg <- h.sp.vec[i] # the species
  gx <- substring(tg, first=1,last=regexpr(tg, pattern = " ")-1) # the genus
  df.ncbi <- ids.for.classif[[i]] # data table from initial NCBI query
  
  if(
    sum(df.ncbi$scientificname == tg) >= 1 && # is there a hit in ncbi for exact species name?
    sum(df.ncbi$rank[which(df.ncbi$scientificname == tg)] == "species") >= 1 # and are the hits at the rank of species?
    ) { 
    # if there was a hit at the rank of species, then save the id and note the query was 'species-level'
    ids <- df.ncbi$uid[which(df.ncbi$scientificname == tg)]
    q.rank <- c(q.rank, rep('species', length(ids)))
    } else {
      # if there was not a hit at the rank of species, query the genus
      ids.for.classif[i] <- get_uid_(gx) # query the genus
      df.ncbi <- ids.for.classif[[i]] # replace old entry in list
      
      ids <- df.ncbi$uid[which(df.ncbi$scientificname == gx & df.ncbi$rank == "genus")] # if hit, save id
      q.rank <- c(q.rank, rep('genus', length(ids))) # note id is genus-level
    }
  
  # save extracted data in vectors
  ncbi.ids <- c(ncbi.ids, ids)
  tg.vec <- c(tg.vec, rep(tg, length(ids)))
  
  if(i == 1) {i.vec <- rep(i,length(ids))} else {i.vec <- c(i.vec, rep(i, length(ids)))}
}

# combine into dataframe
m.df<-data.frame(i.vec, sp.query = tg.vec, q.rank, ncbi.ids)
```

Now, we use the ids to download the classification for each species from NCBI. This also takes awhile.

``` r
ncbi.tax<-classification(unique(ncbi.ids), db = 'ncbi')
save.image("workspace_after_tax_downloads.RData")
```

The returned data are in a list, so we transform them into a dataframe, and select only those taxonomic ranks that are named (e.g. family, order, class, etc.; not 'no rank').

``` r
ncbi.tax <- rbind(ncbi.tax)%>%
  select(name, rank, query)%>%
  filter(rank != "no rank")
```

Then, we re-arrange the data from a long to a wide format, and combine it with the table containing the query info.

``` r
# make sure combos of 'rank' and 'query' do not have duplicate rows, because it causes a problem when 'spreading' the data next
ncbi.tax <- distinct(ncbi.tax)

ncbi.tax<-spread(ncbi.tax, rank, name) #make long df wide

#combine query info in m.df with output from ncbi using ids as keys
ncbi.tax.out <- left_join(m.df, ncbi.tax, by = c("ncbi.ids" = "query"))

#remove non-animal taxonomies (i.e. not metazoans)
ncbi.tax.out <- filter(ncbi.tax.out, kingdom == "Metazoa" | is.na(kingdom))

#remove unneeded objects
rm(df.ncbi, gx, i, i.vec, ids, m.df, q.rank, tg, tg.vec)
```

The proportion of hosts in the database for which there was a hit in the NCBI database was 0.94. However, the downloaded taxonomies are not consistent across all hosts species. For example, when we just look at the classical taxonomic hierarchy (genus, family, order, class, phylum), we see quite a few host species are missing order or class.

``` r
sapply(select(ncbi.tax.out, genus, family, order, class, phylum), function(x) sum(is.na(x)))
```

    ##  genus family  order  class phylum 
    ##      0      5    152    304      2

Instead of ignoring species with missing taxonomic information, we can try to impute it with a quick and dirty algorithm. I created a couple functions to help with this. Essentially, we look to straddling taxonomic ranks to fill in missing ones (e.g. if class is missing, then is there a subclass returned by ncbi?). Sometimes multiple straddling ranks will be available, and in these cases we just return the most common one in the dataset. Essentially, when imputing taxonomic values, we are preferencing commoness. On the one hand, this may reduce the observed taxonomic diversity of hosts for a given parasite, as certain hosts are likely to share the same imputed taxonomies. On the other hand, the observed taxonomic diversity of hosts is probably higher than if species were ignored when they lacked full taxonomic hierarchies.

We need to create a table for use in the algorithm. It records the number of times any given taxonomic name appears in the data, and is used to judge which taxa is the best for filling in missing data.

``` r
common.tax.groups <- select(ncbi.tax.out, class:tribe)%>%
  gather(key = "tax.rank", value = "tax.name", class:tribe)%>%
  group_by(tax.name)%>%
  summarise(count = n())%>%
  na.omit()
```

Then we fill in the missing data, by looping over the taxonomic ranks and the species.

``` r
source("impute_tax_levels.R") # load functions for imputing

ncbi.tax.out2 <- ncbi.tax.out #make a copy to check how much was filled in

tax.ranks <- c('genus', 'family', 'order', 'class', 'phylum')

for(i in seq_along(tax.ranks)){ #the first loop goes across main taxonomic columns in ncbi.tax.out
  
  rank <- tax.ranks[i]
  
  for(j in seq_along(ncbi.tax.out$i.vec)){ #the second loop goes through the rows of main columns
    
    if( is.na(ncbi.tax.out[j, rank])) {
      
      # if a main tax rank is missing, replace it
      replacement_level <- fill_in_tax_level(ncbi.tax.out, ncbi.tax.out[j, 'sp.query'], rank)
      if(length(replacement_level) > 0) { 
        ncbi.tax.out[j, rank] <- replacement_level
        # only replace if replacement available
        # tried to make this conditional superfluous by including it in 'fill_in_tax_level' function
        # but for some reason I couldn't figure out, it sometimes it returns empty character vectors
      }
    }
  }
}
```

Check how much was filled in.

``` r
data.frame(before = sapply(ncbi.tax.out2[, tax.ranks], function(x) sum(is.na(x))), #before
           after = sapply(ncbi.tax.out[, tax.ranks], function(x) sum(is.na(x)))) #after
```

    ##        before after
    ## genus       0     0
    ## family      5     0
    ## order     152     1
    ## class     304     6
    ## phylum      2     2

Quite an improvement. Thus, for species where some taxonomy was returned, then often there was taxonomic ranks outside the classic hierarchy (genus, family, order, class, phylum), making it usually possible to impute missing taxonomic ranks.

Lastly, we clean up the workspace and save the taxonomy table to file.

``` r
rm(common.tax.groups, i, j, ncbi.tax.out2, rank, tax.ranks)

write.table(ncbi.tax.out, file = "../data/ncbi_host_taxonomy.csv", sep = ",", row.names = F) 
```
