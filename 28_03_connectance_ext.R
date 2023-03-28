####3.1 CONNECTANCE ON EXTINCTIONS/ABUNDANCES using a niche-model####

library(ATNr)
library(dplyr)


####3.2 defining explanatory/response variables####  
#explanatory:
#connectance of the food web
#n_sub / n_tot = size of the sub-food web in relation to the total food web

#responses (each to be compared, between a food web of n species and a sub-foodweb of k<n species)
No.ext = NULL #number of extinctions
abundances =NULL #total biomass / Body mass

####3.3 set up parameters and food web####
n_bas <- 8
n_tot <- 32
n_sub = seq(0,1, by=0.125)*(n_tot-n_bas) 


#seeds.fw <- seq(1,length(con)*length(reps), by=1) #seeds for creating the food web, 
#LEAVE OUT FOR NOW

con <- seq(0.1, 0.45, by = 0.05) #connectance values of the food webs (EXPLANATORY)
reps<- 10 #no of replicates per food web
times <- seq(1, 1e8, by = 1e6) #times for integration
biom <- runif(n_tot, 1, 4) #initial biomasses
BM <- runif(n_tot, 2, 3) %>% #body masses of the species
  sort()
BM <- 10^BM



#intialise variables for loop
i = 0
exts = c()

#initialize output vector:
extinctions <- rep(NA, length(con))
biomasses_out <- matrix(data=NA, nrow=length(con), ncol = n_tot) 
  #a matrix containing the biomasses at the final time 
  #for each connectance value (rows) and species (columns)




for (i in 1:length(con)){
  fw <- create_niche_model(S = n_tot, C = con[i])
  
  
  model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
    initialise_default_Unscaled()
  
  #solve ode's
  sol = lsoda_wrapper(times, biom, model)
  exts = sol[nrow(sol), -1]
  exts = sum(exts <= model$ext)
  extinctions[i] = exts
  biomasses_out[i,] = sol[nrow(sol), -1] 
    #takes the biomasses per species and puts them into a matrix
  
}

biomasses_out <- as.data.frame(biomasses_out)
row.names(biomasses_out) <- c(1:8)

plot(biomasses_out$V1)


