library(shinystan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(tidyverse)
library(readxl)
library(here)
library(rethinking)
library(rstan)

set.seed(108)

#function for calculating relative amp efficiency for each taxon
EFFICIENCY <- function(startingDNA, endingReadCount, numberCycles, ...){
  temp <- unlist(
    (endingReadCount/startingDNA)^(1/numberCycles) - 1)
  temp[which(is.infinite(temp))] <- NA
  temp <- temp/(max(temp, na.rm=T))  #normalize to maximum, to create relative index of amp efficiency
  temp[temp<=0] <- 0
  as.numeric(temp)
}

#to fit to a Beta, make strictly non-zero, non-one
unZeroOne <- function(x){
  x[x<1e-5] <- 1e-5  #set zero values to near zero
  x[x> (1 - 1e-5)] <- (1 - 1e-5) #set one values to near one
  x
}

summary.table <- data.frame(
  Reference = c("Port et al. 2016", "Hänfling et al. 2016", "Hänfling et al. 2016", "Olds et al. 2016", "Olds et al. 2016", 
                "Olds et al. 2016", "Olds et al. 2016", "Deiner et al. 2016", "Ford et al. 2016", "Braukmann et al 2019", "Andruszkiewicz et al 2017"),
  GeneRegion = c("12s", "12s", "Cytochrome B", "Cytochrome B", "12s", "16s", "12s", "COI", "16s", "COI", "12s"),
  PrimerName = c("Riaz 12s", "Riaz 12s", "L14841/H15149", "L14841/H15149", "Am12s", "Ac16s", "Ac12s", "Folmer", 
                 "Salmon-F/16s-R", "C_LepFolF/C_LepFolR", "MiFish"),
  TargetGroup = c(rep("Fish", 7), "Eukaryotes", "Fish", "Arthropods", "Fish"),
  SampleSize = c(10,22,22, 6,6,6,6, 33, 5, 374, 10),
  Shape1 = NA,
  Shape1SD = NA,
  Shape2 = NA,
  Shape2SD = NA
)

#################################
dataIn <- read_xlsx(here("Data/PrimerEfficiency.xlsx"))
    #convert proportions to read counts, multiplying by arbitrary number of total reads for run
    dataIn$Ending_reads[dataIn$Author == "Ford"] <- dataIn$Ending_reads[dataIn$Author == "Ford"] * 10^5
#################################


dataOut <- list()
for (i in dataIn$Dataset){
  temp <- dataIn[dataIn$Dataset == i,]
  dataOut[[i]] <- EFFICIENCY(startingDNA = temp[,5], 
                             endingReadCount = temp[,6], 
                             numberCycles = temp[,7])
}


dataForStan <- dataIn %>% 
  unite(AuthorPrimer, c(Author, Primer), sep = "_") %>% 
  mutate(Dataset = match(AuthorPrimer, unique(AuthorPrimer))) %>%  #create dataset index using unique author/primer combinations
  mutate(Taxon = as.numeric(as.factor(Taxon))) %>% 
  mutate(Value = EFFICIENCY(Starting_prop, 
                            Ending_reads,
                            Ncycles)) %>% 
  group_by(Dataset, CommNumber) %>% 
  mutate(Value = Value / max(Value, na.rm = T)) %>%  #scale within PCR reaction; these are therefore relative amp efficiencies
  filter(!is.na(Value)) %>% 
  mutate(Value = unZeroOne(Value)) %>%   #so beta distribution will work  
  select(AuthorPrimer, Dataset, CommNumber, Taxon, Value) %>% 
  mutate(TaxonDataset = match(paste0(AuthorPrimer, Taxon), unique(paste0(AuthorPrimer, Taxon)))) #unique taxon-dataset combinations

#take first Hänfling dataset as example
h <- data.frame(dataOut[[10]], 
           dataOut[[11]],
           dataOut[[12]],
           dataOut[[13]],
           dataOut[[14]],
           dataOut[[15]],
           dataOut[[16]],
           dataOut[[17]],
           dataOut[[18]],
           dataOut[[19]]
)

colnames(h) <- 1:10
rownames(h) <- 1:22

h <- h %>% 
  rownames_to_column("Taxon") %>% 
  gather(key = Community, value = Value, -Taxon) %>% 
  filter(!is.na(Value))

h$Value <- unZeroOne(h$Value) 
h$Taxon <- as.numeric(h$Taxon)
h$Community <- as.numeric(h$Community)

# stanMod <- map2stan(
#   alist(
#     Value ~ dbeta(shape1, shape2),
#     shape1[Community] ~ dlnorm(0,0.3),
#     shape2[Community] ~ dlnorm(0,0.3)
#   ) ,
#   data = h,
#   chains = 5, 
#   iter = 10000,
#   warmup = 500
# )
# stanMod@model

#uses single dataset to take means of amp efficiencies across taxa, w a single shared variance
mymodel <- stan_model(
  file = "Analysis/simple.means.stan")
model.fits <-
  sampling(object = mymodel,
           data = list(
             "Community" = h[,2],
             "Value" = h[,3],
             "Taxon" = h[,1],
             "N" = nrow(h),
             "N_comm" = length(unique(h$Community)),
             "N_tax" = length(unique(h$Taxon)))
  )


launch_shinystan(model.fits)

precis(model.fits, depth = 2)


#The below estimates the mean / variance of amp efficiencies for each taxon, and then
#rolls these up into an estimation of the shape params describing the beta distribution from which 
#they are drawn.  I am trying to propagate error up the hierarchy, rather than reducing variance artificially
#by using sample means to estimate the betas.

BetaWithVariancePropagated <- stan_model(
  file = "Analysis/hierarchicalModelPractice.stan")
Beta.model.fits <-
  sampling(object = BetaWithVariancePropagated,
           data = list(
             #"Community" = h[,2],
             "Value" = h[,3],
             "Taxon" = h[,1],
             "N" = nrow(h),
             #"N_comm" = length(unique(h$Community)),
             "N_tax" = length(unique(h$Taxon)),
             chains = 1)
  )

precis(Beta.model.fits, depth = 2)
launch_shinystan(Beta.model.fits)


#This does the above, but does it across all datasets, and also lets the variance term among replicate amp
#efficiencies for a taxon vary between datasets. The result should be a different set of shape parameters for each 
#dataset

AllDatasetsBetaWithVariancePropagated <- stan_model(
  file = "Analysis/hierarchicalModelPractice_allDatasets.stan")
Beta.model.fits <-
  sampling(object = AllDatasetsBetaWithVariancePropagated,
           data = list(
             "Value" = as.numeric(unlist(dataForStan[,5])),
             "Taxon" = as.numeric(unlist(dataForStan$TaxonDataset)), #taxon-dataset combinations
             "Dataset" = as.numeric(unlist(dataForStan[,2])),
             "N" = nrow(dataForStan),
             "N_tax" = length(unique(dataForStan$TaxonDataset)),  
             "N_datasets" = length(unique(dataForStan$Dataset))
             ),
           chains = 1, iter = 4000
  )

summary(Beta.model.fits)
#launch_shinystan(Beta.model.fits)

Dataset <- unique(as.numeric(unlist(dataForStan[,2]))) #index number of dataset

stan.res <- as.data.frame(summary(Beta.model.fits))[1:(2*length(Dataset)),c(1,3)] #get beta shape params
stan.res <- round(stan.res, 3)

#order to match summary.table
#unique(dataIn[,c("Author", "Primer")])[c(3,8,9,7,6,5,4,2,1,10,11),]
summary.table$Shape1 <- stan.res[Dataset, 1][c(3,8,9,7,6,5,4,2,1,10,11)]
summary.table$Shape2 <- stan.res[(Dataset + length(Dataset)), 1][c(3,8,9,7,6,5,4,2,1,10,11)]
summary.table$Shape1SD <- stan.res[Dataset, 2][c(3,8,9,7,6,5,4,2,1,10,11)]
summary.table$Shape2SD <- stan.res[Dataset + length(Dataset), 2][c(3,8,9,7,6,5,4,2,1,10,11)]


write.csv(summary.table, file = here(paste0("Data/BetaParams_summaryTable_",Sys.Date(),".csv")))

#all of these primer sets do a good job of amplifying the target taxa
par(mfrow = c(6,2))
for (i in 1:nrow(summary.table)){
  hist(rbeta(1000, summary.table$Shape1[i], summary.table$Shape2[i]))
}
  






