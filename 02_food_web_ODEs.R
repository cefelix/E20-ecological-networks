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
ext_thresh <- 0.1**6                    #threshold below which species is considered, extinct


###
####2.2 initialize output arrays####
###

#output arrays: abundances (= densities)
abundances <- array(dim = c(length(con), n_tot))
abundances <- provideDimnames(abundances, sep = "_", base = list( "con", "spec"))

#biomasses
biomasses <- array(dim = c(length(con), n_tot))
biomasses <- provideDimnames(biomasses, sep = "_", base = list( "con", "spec"))

#extinctions
extinctions <- array(dim = c(length(con), n_tot))
extinctions <- provideDimnames(extinctions, sep = "_", base = list( "con", "spec"))

#trophic level
troph.lvl <- array(dim = c(length(con), n_tot))
troph.lvl <- provideDimnames(troph.lvl, sep = "_", base = list( "con", "spec"))

#prey at start (for each species)
prey_on_START <- array(dim = c(length(con), n_tot))
prey_on_START <- provideDimnames(prey_on_START, sep = "_", base = list( "con", "spec"))

#predators at start
consumed_by_START <- array(dim = c(length(con), n_tot))
consumed_by_START <- provideDimnames(consumed_by_START, sep = "_", base = list( "con", "spec"))

#prey at end
prey_on_END <- array(dim = c(length(con), n_tot))
prey_on_END <- provideDimnames(prey_on_END, sep = "_", base = list( "con", "spec"))

#predators at end
consumed_by_END <- array(dim = c(length(con), n_tot))
consumed_by_END <- provideDimnames(consumed_by_END, sep = "_", base = list( "con", "spec"))

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

    BM <- runif(n_tot, 1, 12) %>% 
      sort()
    BM <- (10^BM) #Body masses, spanning 12 orders of magnitude
    
    model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Unscaled()

    model$ext <- ext_thresh
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    
    #calculating output variables
    bioms = sol[nrow(sol),-1]                           #final biomasses
    abuns = bioms / BM                                  #final abundances
    exts = ifelse(bioms <= model$ext, yes = 1, no = 0)  #extinctions
    trolev =  TroLev(fw)                                #trophic levels
    prey.on_start = colSums(fw)                         #trophic links: species is predator
    consumed.by_start = rowSums(fw)                     #trophic links: species is consumed 
    
    #storing them in arrays
    biomasses[i,] = bioms #final biomasses
    abundances[i,] = abuns #abundances
    extinctions[i,] = exts
    troph.lvl[i,] = trolev
    prey_on_START[i,] = prey_links #number of prey species for each species
    consumed_by_START[i,] = pred_links #number of predator species for each species
    
   
    #exclude interaction with extinct species from adjacency matrix
    fw_end <- fw
      fw_end[exts == 1,] <- 0 #all rows which contain an extinct species are set to zero
      fw_end[,exts == 1] <- 0 #all columns which contain an extinct species are set to zero
        sum(colSums(fw) - colSums(fw_end)) #the number of links which are lost due to the extinctions
        sum(rowSums(fw) - rowSums(fw_end)) #equal to above -> everything fine
    
    #calculate number of trophic links at and of simulation    
    prey.on_start = colSums(fw)                       
    consumed.by_start = rowSums(fw) 
    
    #storing them in arrays 
    prey_on_END[i,] = prey.on_start
    consumed_by_END[i,] = consumed.by_start
    
    #count iterations
    print(i)
  }
  




###
####2.4 store output arrays into RDS file####
###
output_96 <- list(abundances, biomasses, extinctions, troph.lvl, 
                  prey_on_START, prey_on_END, consumed_by_START, consumed_by_END)
names(output_96) <- c("abundances", "biomasses", "extinctions", "troph_lvl", 
                      "prey_on_START", "prey_on_END", "consumed_by_START", "consumed_by_END")

saveRDS(output_96, file = "./raw/20220508_96spec_2000c_v01.rds")



###
####2.5 read in RDS file####
###

output_64 <- readRDS(file = "./raw/output_64.rds")

abundances <- output_64$abundances
biomasses <- output_64$biomasses
extinctions <- output_64$extinctions
troph.lvl <- output_64$troph_lvl


###
####2.6 check extinctions####
###
for(i in 1:length(con)) {
  sum(extinctions[1,i,]) %>%
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


