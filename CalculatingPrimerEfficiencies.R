#Calculating primer efficiencies
#Kelly et al. 2019

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

dataIn <- read_xlsx(here("Data/PrimerEfficiency.xlsx"))

#dataIn <-  dataIn %>% 
#    filter(Starting_prop > 0)

#convert proportions to read counts, multiplying by arbitrary number of total reads for run
dataIn$Ending_reads[dataIn$Author == "Ford"] <- dataIn$Ending_reads[dataIn$Author == "Ford"] * 10^5


dataOut <- list()
for (i in dataIn$Dataset){
  temp <- dataIn[dataIn$Dataset == i,]
  dataOut[[i]] <- EFFICIENCY(startingDNA = temp[,5], 
                             endingReadCount = temp[,6], 
                             numberCycles = temp[,7])
}

unique(dataIn[,1:2])
as.data.frame(dataIn[dataIn$Dataset == 9,])

unique(dataIn$Dataset[which(dataIn$Author == "Olds" &
        dataIn$Primer == "L14735_H15149c_cytB")])

#for datasets w multiple replicates of the same primer/community match, take mean of amp efficiencies for each taxon 
ford <- apply(data.frame(dataOut[[1]], dataOut[[2]]), MARGIN = 1, FUN = "mean", na.rm = T)
port <- apply(data.frame(dataOut[[4]], dataOut[[5]]), MARGIN = 1, FUN = "mean", na.rm = T)
Hanfling12s <- apply(data.frame(dataOut[[10]], 
                             dataOut[[11]],
                             dataOut[[12]],
                             dataOut[[13]],
                             dataOut[[14]],
                             dataOut[[15]],
                             dataOut[[16]],
                             dataOut[[17]],
                             dataOut[[18]],
                             dataOut[[19]]
                             ), MARGIN = 1, FUN = "mean", na.rm = T)
HanflingCytB <- apply(data.frame(dataOut[[20]], 
                                dataOut[[21]],
                                dataOut[[22]],
                                dataOut[[23]],
                                dataOut[[24]],
                                dataOut[[25]],
                                dataOut[[26]],
                                dataOut[[27]],
                                dataOut[[28]],
                                dataOut[[29]]), 
                      MARGIN = 1, FUN = "mean", na.rm = T)
dataOut[[30]] <- dataOut[[30]][!is.na(dataOut[[30]])]  #remove NAs

dataOut <- list(
  port,
  Hanfling12s,
  HanflingCytB,
  dataOut[[6]],
  dataOut[[7]],
  dataOut[[8]],
  dataOut[[9]],
  dataOut[[3]],
  ford,
  dataOut[[30]]
)

#make acceptable for beta
dataOut <- lapply(dataOut, unZeroOne)


#use stan to fit beta distr to each 
modelOut <- list()
stanMod <- map2stan(
  alist(
    A ~ dbeta(shape1, shape2),
    shape1 ~ dlnorm(0,0.3),
    shape2 ~ dlnorm(0,0.3)
  ) ,
  data = data.frame(A = dataOut[[9]]),
  chains = 5, 
  iter = 10000,
  warmup = 500
)

  for (i in 1:length(dataOut)) {
    modelOut[[i]] <- map2stan(stanMod, 
                              data = data.frame(A = dataOut[[i]]),
                              chains = 5)
    print(i)
    }
  
i = 1
precis(modelOut[[i]], prob = 0.95)
plot(precis(modelOut[[i]], prob = 0.95))
post <- extract.samples(modelOut[[i]])
dens(post[[1]])
plot(modelOut[[i]])

saveRDS(modelOut, file = "Analysis/stanModelsOut.RDS")

ModelMeans <- data.frame(NA, NA)
ModelSD <- data.frame(NA, NA)
for (i in 1:length(modelOut)){
  ModelMeans[i,] <- precis(modelOut[[i]])@output$Mean
  ModelSD[i,] <- precis(modelOut[[i]])@output$StdDev
}

Model.results <- data.frame(
  unique(dataIn[,1:2]),
  round(ModelMeans, 2),
  round(ModelSD, 2)
)
names(Model.results)[3:6] <- c("Shape1",
                               "Shape2",
                               "Shape1_SD",
                               "Shape2_SD")

summary.table <- data.frame(
  Reference = c("Port et al. 2016", "Hänfling et al. 2016", "Hänfling et al. 2016", "Olds et al. 2016", "Olds et al. 2016", "Olds et al. 2016", "Olds et al. 2016", "Deiner et al. 2016", "Ford et al. 2016", "Braukmann et al 2019"),
  GeneRegion = c("12s", "12s", "Cytochrome B", "Cytochrome B", "12s", "16s", "12s", "COI", "16s", "COI"),
  PrimerName = c("Riaz 12s", "Riaz 12s", "L14841/H15149", "L14841/H15149", "Am12s", "Ac16s", "Ac12s", "Folmer", "Salmon-F/16s-R", "C_LepFolF/C_LepFolR"),
  TargetGroup = c(rep("Fish", 7), "Eukaryotes", "Fish", "Arthropods"),
  SampleSize = c(10,22,22, 6,6,6,6, 33, 5, 374),
  Shape1 = NA,
  Shape1SD = NA,
  Shape2 = NA,
  Shape2SD = NA
)

summary.table[,6:9] <- Model.results[,c(3,5,4,6)]

write.csv(summary.table, "Data/Beta.Stan.Model.Results.csv", row.names = F, quote = F)


