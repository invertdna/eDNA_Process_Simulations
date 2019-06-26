#functions for amplification simulation paper; 2018 - 2019 Kelly Lab

require(vegan)


#function for eDNA sampling, providing stochastic processes for each step in detection/analytical process
ampFunction <- function(B, 
                        ncycles = 35,
                        shedding.eDNA = rlnorm(length(B), 0.5, 0.5), #shedding as function of biomass; output is the eDNA available. Normal-ish distrib
                        a = rbeta(length(B), 1.1, 0.54)  # This is amplification efficiency;parameters taken from 12s fish estimate derived from data in Hänfling et al 2016
) { #where B is a vector of biomass proportions, and ncycles is the number of PCR cycles
  require(vegan)
  
  DNA <- B * shedding.eDNA  # DNA present is biomass times shedding rate
  DNA <- DNA/(sum(DNA))  #rescale, so these are again proportions
  
  epsilon <- rlnorm(n = length(B),  #error term for amplicon generation, adding a bit of stochasticity
                meanlog = 0, 
                sdlog = .05)
  amplicon <- DNA*(((a + 1)^ncycles)*epsilon)  #the number of amplicons present after an `ncycles'-cycle PCR, given amplification efficiency `a` and error term eta
  amplicon <- amplicon/(sum(amplicon))    #rescale to proportions
  
  reads <- rrarefy(
    round(amplicon*10^7, 0),
    round(rbeta(1, 30, 30)*10^6, 0)
  ) #grab random set of 10^5 - 10^6 reads, after assuming 10^7 amplicons were generated, to simulate MiSeq process and pipetting error
  return(reads)
}


#function for calculating various diversity metrics; depends on vegan
diversity.metrics <- function(x) {  #where x is a vector
  list(Richness = specnumber(x),
       Shannon = diversity(x), 
       Simpson = diversity(x, index ="simpson"),
       InvSimpson = diversity(x, index ="invsimpson"),
       Pielou = diversity(x)/log(specnumber(x))
  )
}


#for creating a mix of beta distributions; Ole Shelton 2018
beta_2_mixture <- function(n.species,alpha1,beta1,alpha2,beta2,rho){
  mix <- rbinom(n.species,1,rho)  ## rho is the fraction that comes from the first beta distribution.  The rest comes from the second.
  first  <- rbeta(n.species,alpha1,beta1)
  second <- rbeta(n.species,alpha2,beta2)
  
  all <- first %>% as.data.frame() %>% mutate(val=ifelse(mix==0,second,first)) %>% select(val)
  return(unlist(all))
}


ChangeBiomass <- function(B){
  temp.B <-
    B +  #biomass vector in
    (sample(c(-1,1), size = length(B), replace = T) * # probabilistic increase or decrease of...
       rbeta(length(B), 5, 5) * #a percentage drawn from a normal-ish beta distribution of mean 0.5
       B)
  temp.B/sum(temp.B)  #renormalize to proportions
}


ChangeBiomassDirichlet <- function(B){
  as.vector(rdirichlet(n = 1, alpha = rep(5, times = length(B)))) #random biomass with each time step, rather than being an autoregressive process
}

#Permutation function; by column and row
PERMFUN <- function(z)    {
  z <- as.matrix(z)
  as.data.frame(matrix(sample(z),nrow=nrow(z)))
} #permute function


#get mode of a distribution
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#make elements of a list into a tidy df
MAKETIDY <- function(LIST, index = 1) {
  require(tidyverse)
  LIST[[index]] %>% 
    as.data.frame %>% 
    rownames_to_column("Taxon") %>% 
    gather(key = "Timepoint", value = !!names(LIST)[[index]], -Taxon)
}

#get equation of ggplot lm ; modified from https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph 
lm_eqn <- function(df){
  y <- df[,1]
  x <- df[,2]
  m <- lm(y ~ x)
  eq <- substitute(~~r^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  
  as.character(as.expression(eq))                 
}

#custom log function for dealing w zeros
mylog <- function(x){
  sapply(x, 
         FUN = function(z) if (z == 0) 0 else log(z, base = 2) + 1
  )
}

#function to calculate eDNA index
eDNAINDEX <- function(x) { #where x is a dataframe with taxa/OTUs/etc in rows, and samples in columns
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}
#df <- as.data.frame(matrix(sample(1:10, 100, replace = T), ncol = 10, nrow = 10))
#eDNAINDEX(df)

#get max of a row
rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}


#a function to standardize a dataframe of eDNA counts (rows = OTUs or similar, columns = samples) 
#by all of the methods available in vegan's decostand function; 
#output is a list of standardized dfs
ALLSTAND <- function(df) {
  require(vegan)
  res.list <- list()
  methodvec <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")
  marginvec <- c(2,1,2,2,1,2,1,NA,2,2,NA) #to override default margins in vegan
  
  for (i in 1:length(methodvec)){
    res.list[[i]] <- decostand(df, 
                               method = methodvec[i],
                               MARGIN = marginvec[i])
  }
  names(res.list) <- methodvec
  return(res.list) 
}


#Calculate Spearman's Rho in dplyr / tidyverse notation; useful for applying to subsets of data
TIDYSPEARMAN <- function(index.results = index.results, 
                         metric = "eDNA_reads"){
  require(vegan)
  metricquo <- enquo(metric)  #see https://dplyr.tidyverse.org/articles/programming.html
  metricvar <- as.symbol(metric)
  index.results %>% 
    mutate(Taxon = as.character(Taxon)) %>% 
    group_by(Taxon) %>% 
    do(Rho = cor(.[,"Biomass"],
                 .[,!!metricquo], 
                 method = "spearman")) %>% 
    mutate(Rho = as.vector(Rho)) %>% 
    mutate(Statistic = !!metricquo)
}


# Multiple plot function; from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#function for calculating relative amp efficiency for each taxon
EFFICIENCY <- function(startingDNA, endingReadCount, numberCycles){
  temp <- (endingReadCount/startingDNA)^(1/numberCycles) - 1
  temp <- temp/max(temp, na.rm=T)  #normalize to maximum, to create relative index of amp efficiency
  temp[temp<=0] <- 0 
  temp
}

