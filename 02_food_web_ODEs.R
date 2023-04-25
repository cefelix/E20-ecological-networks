####0.0 libraries####
#test

library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


####2 CONNECTANCE ON EXTINCTIONS / ABUNDANCES / SHANNON INDEX, using a niche-model####
###
###
####2.1 set up parameters and food web####
#number of basal species depends on the food web!
n_tot <- 128              #no. of all species

perc_sub <- 0.25          #share of the species present in sub-food webs (float between 0 and 1)
n_sub <- perc_sub*n_tot   #no. of species in a sub-foodweb

con <- seq(0.025, 0.05, by = 0.00125)   #connectance of the food webs 
                                        #check out https://www.pnas.org/doi/epdf/10.1073/pnas.192407699
                                        #mininmal connectance 0.026, maximal 0.315 (from 17 empirical food webs, fig.1)

reps<- 5                                #no of replicates per food web
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
abundance_array <- array(dim = c(reps,length(con), n_tot))
abundance_array <- provideDimnames(abundance_array, sep = "_", base = list("rep", "con", "spec"))

#biomasses
biomass_array <- array(dim = c(reps,length(con), n_tot))
biomass_array <- provideDimnames(biomass_array, sep = "_", base = list("rep", "con", "spec"))

#extinctions
extinction_array <- array(dim = c(reps,length(con), n_tot))
extinction_array <- provideDimnames(extinction_array, sep = "_", base = list("rep", "con", "spec"))

#trophic level
troph.lvl_array <- array(dim = c(reps,length(con), n_tot))
troph.lvl_array <- provideDimnames(troph.lvl_array, sep = "_", base = list("rep", "con", "spec"))



#loop counters (i set them up against the alphabet because that's cooler)
j=1
i=1


####2.3 compute 1e10 food webs####
for (j in 1:reps) {
  BM <- runif(n_tot, 1, 12) %>% #body masses of the species
    sort()
  BM <- (10^BM)
  i=0
  
  for (i in 1:length(con)){
    fw <- create_niche_model(S = n_tot, C = con[i])
    n_bas <- sum(colSums(fw) == 0)
    
    model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Unscaled()
    model$ext <- ext_thresh
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    
    #storing abundances in the abundance array
    abuns = sol[nrow(sol),-1] / BM
    abundance_array[j,i,] = abuns
    
    #final biomasses
    bioms = sol[nrow(sol),-1]
    biomass_array[j,i,] = bioms
    
    #extinctions
    exts = sol[nrow(sol), -1]
    exts = ifelse(exts <= model$ext, yes = 1, no = 0) #accordingly, 1 means species is extinct!
    extinction_array[j,i,] = exts
    
    #trophic levels 
    troph.lvl_array[j,i,] = TroLev(fw)
    
  }
  
  print(j)
}


###
####2.4 store output arrays into RDS file####
###

#DANGER: this stores the results from the ODEs above 
saveRDS(abundance_array, file = "./raw/abundance_conHRs.rds")
saveRDS(biomass_array, file = "./raw/biomass_conHR.rds")
saveRDS(extinction_array, file = "./raw/extinction_conHR.rds")
saveRDS(troph.lvl_array, file = "./raw/trophic_lvl_conHR.rds")


###
####2.5 read in RDS files####
###

#read RDS files:
abundance_array <- readRDS(file = "./raw/abundance.rds")
biomass_array <- readRDS(file = "./raw/biomass.rds")
extinction_array <-  readRDS(file = "./raw/extinction.rds")
troph.lvl_array <- readRDS(file = "./raw/trophic_lvl.rds")

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


