####0.0 libraries####

library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


####2 CONNECTANCE ON EXTINCTIONS / ABUNDANCES / SHANNON INDEX, using a niche-model####
###
###
####2.1 set up parameters and food web####
#number of basal species depends on the food web!
n_tot <- 96            #no. of all species

#perc_sub <- seq(1/16, 1/4, by=1/32)  
perc_sub <- 0.25          #share of the species present in sub-food webs (float between 0 and 1)
n_sub <- perc_sub*n_tot   #no. of species in a sub-foodweb

set.seed(666)
#con <- seq(0.025, 0.05, by = 0.025)   #connectance of the food webs 
                                        #check out https://www.pnas.org/doi/epdf/10.1073/pnas.192407699
                                        #mininmal connectance 0.026, maximal 0.315 (from 17 empirical food webs, fig.1)
con <- runif(2000, min=0.05, max=0.35) %>% sort()


times <- seq(1, 1e12, by = 1e9)         #time for integration of dynamics
biom <- runif(n_tot, 1, 4)              #initial biomasses
BM <- runif(n_tot, 1, 12) %>%           #Body masses 
  sort()                                #realistic spanning range in soil food webs (Potapov 2021):
BM <- (10^BM)                           #https://pubmed.ncbi.nlm.nih.gov/34086977/
ext_thresh <- 0.1**6                    #threshold below which species is considered, extinct


###
####2.2 initialize output arrays####
###

#output arrays: abundances (= densities)
abundance_array <- array(dim = c(length(con), n_tot))
abundance_array <- provideDimnames(abundance_array, sep = "_", base = list( "con", "spec"))

#biomasses
biomass_array <- array(dim = c(length(con), n_tot))
biomass_array <- provideDimnames(biomass_array, sep = "_", base = list( "con", "spec"))

#extinctions
extinction_array <- array(dim = c(length(con), n_tot))
extinction_array <- provideDimnames(extinction_array, sep = "_", base = list( "con", "spec"))

#trophic level
troph.lvl_array <- array(dim = c(length(con), n_tot))
troph.lvl_array <- provideDimnames(troph.lvl_array, sep = "_", base = list( "con", "spec"))

#prey at start (for each species)
prey_array <- array(dim = c(length(con), n_tot))
prey_array <- provideDimnames(prey_array, sep = "_", base = list( "con", "spec"))

#predators at start
predators_array <- array(dim = c(length(con), n_tot))
predators_array <- provideDimnames(predators_array, sep = "_", base = list( "con", "spec"))

#prey at end
prey_array.end <- array(dim = c(length(con), n_tot))
prey_array.end <- provideDimnames(prey_array.end, sep = "_", base = list( "con", "spec"))

#predators at end
predators_array.end <- array(dim = c(length(con), n_tot))
predators_array.end <- provideDimnames(predators_array.end, sep = "_", base = list( "con", "spec"))

#loop counters (i set them up against the alphabet because that's cooler)
j=1
i=1


####2.3 compute 1e10 food webs####
  for (i in 1:length(con)){
    
    fw <- create_niche_model(S = n_tot, C = con[i])
    n_bas <- sum(colSums(fw) == 0)
    
    while (n_bas > 12 | n_bas < 4) { #if fw has less than 4 or more than 16 basal species, repeat
      fw <- create_niche_model(S = n_tot, C = con[i])
      n_bas <- sum(colSums(fw) == 0)
    }
    
    BM <- runif(n_tot, 1, 12) %>% #body masses of the species
      sort()
    BM <- (10^BM)
    
    model <- create_model_Scaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Scaled()
    model$ext <- ext_thresh
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    
    #storing abundances in the abundance array
    abuns = sol[nrow(sol),-1] / BM
    abundance_array[i,] = abuns
    
    #final biomasses
    bioms = sol[nrow(sol),-1]
    biomass_array[i,] = bioms
    
    #extinctions
    exts = sol[nrow(sol), -1]
    exts = ifelse(exts <= model$ext, yes = 1, no = 0) #accordingly, 1 means species is extinct!
    extinction_array[i,] = exts
    
    #trophic levels 
    troph.lvl_array[i,] = TroLev(fw)
    
    #preys
    prey_array[i,] = colSums(fw)
    
    #predators
    predators_array[i,] = rowSums(fw)
    
   
    #to do: shrink food web prior to rowsum/colsum calculation:
    fw_end <- fw
    fw_end[exts == 1,] <- 0 #all rows which contain an extinct species are set to zero
    fw_end[,exts == 1] <- 0 #all columns which contain an extinct species are set to zero
      sum(colSums(fw) - colSums(fw_end)) #the number of links which are lost due to the extinctions
      sum(rowSums(fw) - rowSums(fw_end)) #equal to above -> everything fine
    
    #preys at end
    prey_array.end[i,] = colSums(fw_end)
    predators_array.end[i,] = rowSums(fw_end)
    print(i)
  }
  




###
####2.4 store output arrays into RDS file####
###
output_96 <- list(abundance_array, biomass_array, extinction_array, troph.lvl_array, 
                  prey_array, prey_array.end, predators_array, predators_array.end)
names(output_96) <- c("abundances", "biomasses", "extinctions", "troph_lvl", 
                      "feed_on_START", "feed_on_END", "consumed_by_START", "consumed_by_END")
saveRDS(output_96, file = "./raw/20220508_96spec_2000c_v01.rds")


###
####2.5 read in RDS file####
###

output_64 <- readRDS(file = "./raw/output_64.rds")

abundance_array <- output_64$abundances
biomass_array <- output_64$biomasses
extinction_array <- output_64$extinctions
troph.lvl_array <- output_64$troph_lvl


###
####2.6 check extinctions####
###
for(i in 1:length(con)) {
  sum(extinction_array[1,i,]) %>%
    print()
}


####notes from lecture####

#Neo Martinez
#model doesnt converge: change parameters until it converges


#library(igraph) 
#plot() arguments:
#edges: links
#vertex: nodes
#https://kateto.net/netscix2016.html


###
#purpose of protocol:
#base for grading
#include information i might need in the future to handle similar problems
#how does question translate into an equation


