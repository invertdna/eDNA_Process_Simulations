#Generate Plots for eDNA Simulation Paper -- Kelly et al. 2019

library(vegan, quietly = T)
library(here)
library(broom)
library(tidyverse, quietly = T)
library(patchwork)
library(RColorBrewer)
library(GGally, quietly = T)

source(here("Functions/FunctionsForSimulationPaper.R"))

load(here("simulations.workspace.Rdata"))  #load in data from main analysis Rmd file

#note to future self, and anyone else reading this: I called Case B `left-skewed` in the variable names,
#because the mode is on the left.  But Ole set me straight: in math terms, this is called `right-skewed`, 
#because the median is to the right of the mode. So I changed the figure captions, 
#but left the variable names as `left.skewed`, for ease of use w the main data script.


######Figure 1: Biomass distributions
      pdf("Figures/BiomassDistributions.pdf", width = 7, height = 4 )
        data.frame(B.uniform, 
                 as.vector(B.lessVariable),
                 as.vector(B.moreVariable)) %>% 
        gather(key = "Dataset", value = "value") %>% 
        ggplot(aes(x = value, fill = Dataset, color = Dataset)) + 
        geom_density(position = "identity", alpha = .5, bw = .0004, trim = TRUE) +
        xlim(c(0, .005)) +
        ylim(c(0,1000)) +
        xlab("Proportion of Biomass")  +
        ylab("Number of Taxa") +
        scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(2,4,6)],
                          labels = c("Moderately Variable \n (Gamma = 5)", "More Variable \n(Gamma = 1)", "Uniform"),
                          name = "Biomass Distribution") +
        scale_color_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(2,4,6)],
                           labels = c("Moderately Variable \n (Gamma = 5)", "More Variable \n(Gamma = 1)", "Uniform"),
                           name = "Biomass Distribution") +
        theme_bw() +
        theme(legend.position = c(.75, 0.35)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      dev.off()

      
      
######Figure 2: amplification efficiency distributions
      CaseB.left.skewed.plot <- ggplot(data = as.data.frame(CaseB.left.skewed), 
                                   aes(CaseB.left.skewed)) + 
                                    geom_density(aes(y = ..scaled..), bw = 0.05, 
                                                 fill = RColorBrewer::brewer.pal(6, "PuBu")[c(4)], alpha = .5,
                                                 color = RColorBrewer::brewer.pal(6, "PuBu")[c(4)]) + 
                                    annotate("text", x = 0.9, y = 0.9, label = "B", size = 10) +
                                    xlab("") + xlim(c(0,1)) +
                                    ylab("Probability") + theme_bw() +
                                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        CaseA.approx.normal.plot <- ggplot(data = as.data.frame(CaseA.approx.normal), 
                                    aes(CaseA.approx.normal)) + 
                                    geom_density(aes(y = ..scaled..), bw = 0.05, 
                                                 fill = RColorBrewer::brewer.pal(6, "PuBu")[c(4)], alpha = .5,
                                                 color = RColorBrewer::brewer.pal(6, "PuBu")[c(4)]) + 
                                    annotate("text", x = 0.9, y = 0.9, label = "A", size = 10) +
                                    xlab("") + xlim(c(0,1)) +
                                    ylab("Probability") + theme_bw() +
                                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        CaseC.riaz12s.plot <- ggplot(data = as.data.frame(CaseC.riaz12s), 
                               aes(CaseC.riaz12s)) + 
                              geom_density(aes(y = ..scaled..), 
                                           fill = RColorBrewer::brewer.pal(6, "PuBu")[c(4)], bw = 0.1, alpha = .5,
                                           color = RColorBrewer::brewer.pal(6, "PuBu")[c(4)]) + 
                              annotate("text", x = 0.9, y = 0.9, label = "C", size = 10) +
                              xlab("Amplification Efficiency") + xlim(c(0,1)) +
                              ylab("Probability") + theme_bw() +
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        

        pdf("Figures/AmplificationDistributions.pdf")
        multiplot(CaseA.approx.normal.plot, CaseB.left.skewed.plot, CaseC.riaz12s.plot)
        dev.off()


######Figure 3: Effect of PCR cycles
        
              cycles.Richness.plot <- Tidy.comm.sim.results %>% 
                filter(Index %in% c("Richness")) %>% 
                filter(Origin == "eDNA") %>% 
                ggplot(aes(y = Index.value, x = PCR.cycles, color = amp.Distribution)) +
                geom_point(alpha = .1) +
                stat_smooth(se = F, span = .55, method = "loess", alpha = .5) +
                annotate("text", x = 47, y = 1040, label = "A", size = 10) +
                xlab("PCR Cycles") + ylab("Richness") +
                scale_color_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(3:6)],
                                   guide = F) +
                theme_bw()+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              
              
              cycles.Shannon.plot <- Tidy.comm.sim.results %>% 
                filter(Index %in% c("Shannon")) %>% 
                filter(Origin == "eDNA") %>% 
                ggplot(aes(y = Index.value, x = PCR.cycles, color = amp.Distribution)) +
                geom_point(alpha = .1) +
                stat_smooth(se = F, span = .55, method = "loess", alpha = .5) +
                annotate("text", x = 47, y = 6.5, label = "B", size = 10) +
                xlab("PCR Cycles") + ylab("Shannon Index") +
                scale_color_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(3:6)],
                                   labels = c("Case A : \n  Symmetrical", "Case B : \n  Right-Skewed", "Case C : \n  Limited-Target"),
                                   name = "Amplification \nEfficiency") +
                theme_bw()+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(r = 1, unit = "cm"))
              
        pdf("Figures/SummaryStats_by_pcrCycles_and_ampDistribution.pdf", width = 12, height = 6)
          cycles.Richness.plot + cycles.Shannon.plot  #using patchwork to do multiplot
        dev.off()  


######Figure 4: Effect of Amp Efficiency

        ampRichness.plot <-   Tidy.eff.sim.results %>% 
          filter(Origin == "eDNA", Index == "Richness") %>% 
          ggplot(aes(x= amp.Distribution, y = Index.value, fill = BiomassDistribution))+
          geom_boxplot(outlier.size = 0) +
          #facet_grid(amp.Distribution~., scales = "free_y") +
          scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(4:6)], name = "Biomass\nDistribution", guide = FALSE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_x_discrete(labels = c("Case A : \n  Symmetrical", "Case B : \n  Right-Skewed", "Case C : \n  Limited-Target")) +
          ylab("Richness") + xlab("") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        ampShannon.plot <-     Tidy.eff.sim.results %>% 
          filter(Origin == "eDNA", Index == "Shannon") %>% 
          ggplot(aes(x= amp.Distribution, y = Index.value, fill = BiomassDistribution))+
          geom_boxplot(outlier.size = 0) +
          #facet_grid(amp.Distribution~., scales = "free_y") +
          scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(4:6)], 
                            name = "Biomass\nDistribution",
                            labels = c("Moderately Variable \n (Gamma = 5)", "More Variable \n(Gamma = 1)", "Uniform")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_x_discrete(labels = c("Case A : \n  Symmetrical", "Case B : \n  Right-Skewed", "Case C : \n  Limited-Target")) +
          ylab("Shannon Index") + xlab("") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        
        pdf("Figures/EffectAmpEfficiency_Richness_Shannon.pdf", width = 10, height = 6)
        ampRichness.plot + ampShannon.plot   #using patchwork to create multiplot
        dev.off()  




        
######Figure 5: Probability of Detection
        
        pdf("Figures/LikelihoodDetection_by_AmpEfficiency_Biomass.pdf")
        Detection %>% 
          ggplot(aes(y = counts, x = ampEff, color = biomass)) +
          geom_point(size = .5) +
          stat_smooth(method="glm", method.args=list(family="binomial"), se=TRUE) +
          scale_color_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[c(2,4,6)],
                             labels = c("Uniform", "Less Variable", "More Variable"),
                             name = "Biomass Distribution") +
          ylab("Likelihood of Detection") + xlab("Amplification Efficiency") +
          theme_bw() +
          theme(legend.position = c(.75, 0.25)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        dev.off()


        
######Figure 7: Index performance
        
 
    pdf(here("Figures/IndexPerformance.pdf"), height = 10, width = 7)
      index.results %>% 
        group_by(Taxon) %>% 
        summarise(read.mean = mean(eDNA.total.read.number, na.rm = T)) %>% 
        left_join(tidy.cor) %>% 
        left_join(within.timepoint.Rho) %>% 
        mutate(log.read.mean.quartile = cut_number(log(read.mean), n =4, labels = F)) %>%
        mutate(log.read.mean.quartile = as.factor(log.read.mean.quartile)) %>% 
        filter(!is.na(read.mean)) %>% 
        mutate(Statistic = factor(Statistic, levels = StatOrder)) %>% 
        group_by(Statistic) %>% 
        mutate(medianRho = median(Rho), modeRho = getmode(Rho)) %>% 
        filter(!Statistic %in% c("total", "chi.square", "pa", "range", "standardize", "max", "eDNA_root_index")) %>%
        ggplot(aes(Rho)) + 
        geom_histogram(aes(fill =log.read.mean.quartile), color = NA, binwidth = 0.02) + 
        geom_vline(aes(group = Statistic, xintercept = medianSingleTimeRho + 0.015), color = "grey") + 
        geom_vline(aes(group = Statistic, xintercept = medianRho + 0.015), color = "red") + 
        geom_vline(aes(group = Statistic, xintercept = modeRho + 0.015), color = "darkred") + 
        facet_grid(Statistic ~ ., scales = "free_y") +
        scale_fill_manual(values = RColorBrewer::brewer.pal(6, "PuBu")[3:6]) +
        ylab("Count") + 
        xlab("Spearman's Rho Correlation with Biomass") +
        guides(fill=guide_legend(title="Amplicon\nFrequency\nQuartile")) +
        theme(legend.title = element_text(size=10)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlim(c(0,1)) #note this cuts off half of the null distribution, but that is nonsense anyway, since it's a negative correlation.
      dev.off()

#Figure 8: Rho by Amplification Efficiency

pdf(here("Figures/Rho_by_AmpEfficiency.pdf"), height = 7, width = 5)
      index.results %>% 
        group_by(Taxon) %>% 
        summarise(ampEfficiency = mean(ampEfficiency), Biomass = mean(Biomass)) %>% 
        mutate(`Relative Mean Biomass` = Biomass / min(Biomass)) %>% 
        left_join(tidy.cor) %>% 
        mutate(Statistic = factor(Statistic, levels = StatOrder)) %>% 
        filter(Statistic == "eDNA Index") %>% 
        ggplot(aes(x = ampEfficiency, y = Rho, color = `Relative Mean Biomass`)) +
        geom_point() +
        xlab("Amplification Efficiency") +
        scale_color_continuous(name = "Relative Mean Biomass") +
        theme_bw() +
        theme(legend.position = c(.75, 0.25)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()      
      
      
      
#########################################
#####SUPPLEMENTAL PLOTS
#########################################      

####Port et al. amplification efficiencies
    
        port<-read.csv("Data/Port_et_al_2016_ampEfficiency.csv")
        
        DNA1 <- port[,3] #starting concentrations of DNA for two different mixes of 10 species (community 1 has equal proportions, community 2 varies by a factor of 20)
        DNA2 <- port[,5]
        amplicon1 <- port[,4]  #recovered mean amplicon counts for each species in each community
        amplicon2 <- port[,6]
        ncycles <- 40  #number of PCR cycles, as provided by the paper's methods
        
        com1 <- (amplicon1/DNA1)^(1/ncycles) - 1   #calculating amp efficiency for each taxon
        com2 <- (amplicon2/DNA2)^(1/ncycles) - 1
        
            #scale to max of each
            com1 <- com1/max(com1)
            com2 <- com2/max(com2)
        
        com1[com1 < 0] <- 0      #adjusting small negative number to zero, because negative efficiencies are impossible and this is a result of sampling error
        
        
        pdf(here("Figures/Suppl_Port_et_al.pdf"))
        data.frame(port$Taxon, com1, com2) %>% 
          ggplot(aes(x = com2, y =com1)) +
          geom_point() +
          geom_smooth(method = "lm") +
          xlab("Community 2, Relative Amp. Efficiency") +
          ylab("Community 1, Relative Amp. Efficiency") +
          ylim(c(0,1)) +
          theme_bw() +
          theme_bw()+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        dev.off()


####Hänfling et al. amplification efficiencies
        
        
        hanfling <- read.csv("Data/Hanfling_12s_data.csv")
        
        ##starting concentrations of DNA for different communities
        
        DNA <- hanfling[,seq(3,21,2)]
        names(DNA) <- paste0("DNA_MC", seq(1, 10, 1))
        DNA[DNA==0]<-NA  #DNA not in the original mock community cannot be amplified
        
        amplicon <- hanfling[,seq(2,20,2)]
        names(amplicon) <- paste0("reads_MC", seq(1, 10, 1))
        
        ncycles <- 40  #number of PCR cycles, as provided by the paper's supplemental methods
        
        eff.results <- as.data.frame(matrix(NA, nrow=nrow(DNA), ncol=ncol(DNA)))
        for (i in 1:ncol(DNA)){
          eff.results[,i] <- EFFICIENCY(DNA[,i], 
                                        amplicon[,i],
                                        ncycles)
        }
        colnames(eff.results) <- paste0("MC", seq(1,10,1))
        
        pdf(here("Figures/Suppl_Hanfling_12s.pdf"))
        ggpairs(eff.results)  #within-taxon efficiencies are highly consistent across mock communities w different mixes of species.
        dev.off()
        
        
        #import data
        hanfling <- read.csv("Data/Hanfling_cytB_data.csv")
        
        ##starting concentrations of DNA for different communities
        
        DNA <- hanfling[,seq(3,21,2)]
        names(DNA) <- paste0("DNA_MC", seq(1, 10, 1))
        DNA[DNA==0]<-NA  #DNA not in the original mock community cannot be amplified
        
        amplicon <- hanfling[,seq(2,20,2)]
        names(amplicon) <- paste0("reads_MC", seq(1, 10, 1))
        
        ncycles <- 40  #number of PCR cycles, as provided by the paper's supplemental methods
        
        
        eff.results <- as.data.frame(matrix(NA, nrow=nrow(DNA), ncol=ncol(DNA)))
        for (i in 1:ncol(DNA)){
          eff.results[,i] <- EFFICIENCY(DNA[,i], 
                                        amplicon[,i],
                                        ncycles)
        }          

        colnames(eff.results) <- paste0("MC", seq(1,10,1))
        
        pdf(here("Figures/Suppl_Hanfling_CytB.pdf"))
        ggpairs(eff.results)  #within-taxon efficiencies are highly consistent across mock communities w different mixes of species.
        dev.off()



#######Accumulation Curve by PCR cycle number
      
      customColors <- colorRampPalette(brewer.pal(9, "BuPu")[3:9])(10) #make a 10-color palette
    
pdf("Figures/AccumulationCurve_Ncycles.pdf")
        accumm.sim.curves %>% 
        ggplot(aes(x = as.numeric(SampleSize), y = as.numeric(Ntaxa), color = as.ordered(Ncycles))) +
        geom_line(size = 1.3) +
        xlim(c(0, 25000)) +
        xlab("Sample Size (Reads") + ylab("Species Detected") + 
        scale_color_manual(values = customColors, name = "PCR Cycles") +
        theme_bw() +
        theme(legend.position = c(.75, .1)) +
        guides(color=guide_legend(ncol=5))
dev.off()      
      

#######Figure Multilocus Index by Biomass
      
      
      pdf(here("Figures/Supp_MultilocusIndexPerformance.pdf"), height = 7, width = 10)
      joint.tidy %>% 
        TIDYSPEARMAN(index.results = ., 
                    metric = "eDNA_index") %>% 
        group_by(Statistic) %>% 
        mutate(medianRho = median(Rho), modeRho = getmode(Rho)) %>% 
        ggplot(aes(Rho)) + 
        geom_histogram(binwidth = 0.05, 
                       fill = RColorBrewer::brewer.pal(6, "PuBu")[6], alpha = 0.9) +
        geom_vline(aes(xintercept = medianRho ), color = "red", size = 1.2) + 
        geom_vline(aes(xintercept = modeRho ), color = "darkred", size = 1.2) + 
        ylab("Count") + 
        xlab("Spearman's Rho Correlation with Biomass") +
        theme_bw() +
        theme(legend.title = element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlim(c(0,1)) #note this cuts off half of the null distribution, but that is nonsense anyway, since it's a negative correlation.
      dev.off()
      
