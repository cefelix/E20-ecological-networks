####0.0 libraries####

library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


####2 CONNECTANCE ON EXTINCTIONS / ABUNDANCES / SHANNON INDEX, using a niche-model####
###
###
####2.1 set up parameters and food web####

n_tot <- 96                             #no. of all species
n_sub <- 0.25*n_sub                     #no. of species in a sub-foodweb

set.seed(666)                           
con <- runif(2000, min=0.05, max=0.35) 
con <- sort(con)                        #generating connectance as a continuous variable

times <- seq(1, 1e12, by = 1e9)         #time for integration of dynamics
biom <- runif(n_tot, 1, 4)              #initial biomasses
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




####2.3 compute food web dynamics####
i=1
for (i in 1:length(con)){
  
  BM <- runif(n_tot, 1, 12) %>% #body masses of the species
    sort()
  BM <- (10^BM)
  
  fw <- create_niche_model(S = n_tot, C = con[i])
  n_bas <- sum(colSums(fw) == 0)
  
  while (n_bas > 16 | n_bas < 4) { #if fw has less than 4 or more than 16 basal species, repeat
    fw <- create_niche_model(S = n_tot, C = con[i])
    n_bas <- sum(colSums(fw) == 0)
  }
  
  
  model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>% #creating a niche-model
    initialise_default_Unscaled()
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
  
  #preys for each species at end
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
saveRDS(output_96, file = "./raw/20220505_96spec_15cons_v01.rds")




###
####1 load raw RDS files, initialize parameters####
###

#1.1.1 read RDS list
d.list <- readRDS("./raw/20220505_96spec_15cons_v01.rds")


abundances <- d.list$abundances
biomasses <- d.list$biomasses
extinctions <- d.list$extinctions
troph.lvl <- d.list$troph_lvl

food_START <- d.list$feed_on_START
food_END <- d.list$feed_on_END
consumed_bySTART <- d.list$consumed_by_START
consumed_byEND <- d.list$consumed_by_END


#1.2 initialize selection variables


#first initialize model parameters: go to  ODE-script
slct_bi <- c((n_tot-n_sub+1):n_tot)       #selects biggest consumer species (n_sub=32) 


###
####2 response variable matrices####
### 

#2.1 initializing output arrays
n.bas_mat <- array(NA, dim = c( length(con) ))


extinctions_mat <- apply(extinctions, MARGIN = c(1), FUN = sum) %>% #sum of extinctions for each connectance/replicate combination
  as.matrix()
extinctions.small_mat <- array(NA, dim = c( length(con) ))
extinctions.BAS_mat <- array(NA, dim = c( length(con) ))
extinctions.big_mat <- array(NA, dim = c( length(con) ))




shannon_mat <-  apply(abundances, MARGIN = c(1), FUN = diversity) %>%
  as.matrix() #whole food web
shannon.small_mat <- array(NA, dim = c( length(con) ))
shannon.BAS_mat <- array(NA, dim = c( length(con) ))
shannon.big_mat <- array(NA, dim = c( length(con) ))


rand_select <- seq(from = 1.5*n_sub, to = 4, by=-4) #the sample size for the random selections
shannon.rand_mat        <- array(NA, dim = c( length(con), length(rand_select) )) #randomly selected 
extinctions.rand_mat    <- array(NA, dim = c( length(con), length(rand_select) ))

shannon.rand_matTL      <- array(NA, dim = c( length(con), length(rand_select) )) #randomly selected, but from different trophic levels
extinctions.rand_matTL  <- array(NA, dim = c( length(con), length(rand_select) ))


#2.2 calculating extinctions and shannon indices

slct_bi <- c((n_tot-n_sub+1):n_tot)       #selects biggest consumer species (n_sub=32) 
slct_sm <- c()


for(j in 1:length(con)) {
  
  bas <- n_tot - sum(troph.lvl[j,] > 1)             # the amount of basal species
  consumers <- c((bas+1):n_tot)                       # a vector of all consumer species
  n.bas_mat[j] <- bas                               # the number of basal species in the system
  
  slct_sm <- consumers[1:n_sub]                       # the smallest consumer species
  slct_bi <- c((n_tot-n_sub+1):n_tot)                 # the biggest consumer species
  
  
  #extinctions in 24 species
  extinctions.small_mat[j] <- sum(extinctions[j, slct_sm])
  # selects extinctions in the smallest consumer species
  extinctions.BAS_mat[j] <- sum(extinctions[j, c(1:bas)])
  # selects extinctions in the basal species
  extinctions.big_mat[j] <- sum(extinctions[j, slct_bi])
  # extinctions in the biggest consumer species
  
  #diversity in 24 species
  shannon.small_mat[j] <- diversity(abundances[j, slct_sm])
  # calculates shannon index in the smallest consumer species
  shannon.BAS_mat[j]   <- diversity(abundances[j, c(1:bas)])
  # calculates shannon index in the basal species
  shannon.big_mat[j] <- diversity(abundances[j, slct_bi])
  # shannon index in the biggest consumer species
  
  
  #random selection of species
  for (k in 1:length(rand_select)) {
    slct_rand <- sample(consumers, size = rand_select[k]) %>%       # a random selection of consumer species
      sort()
    
    extinctions.rand_mat[j,k] <- sum(extinctions[j, slct_rand])     # sum of extinctions at con=j for k selections of random species
    shannon.rand_mat[j,k] <- diversity((abundances[j, slct_rand]))  #shannon index at con=j
  }
  
  #random selection of species, but from different trophic levels
  for (k in 1:length(rand_select)) {
    slct_rand <- sample(x = consumers[troph.lvl[j,consumers] < quantile(troph.lvl[j,consumers], 0.25)], #lowest 25% of Trophic levels (consumers)
                        size = rand_select[k]/4) 
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[(troph.lvl[j,consumers] > quantile(troph.lvl[j,consumers], 0.25)) &
                                      (troph.lvl[j,consumers] < quantile(troph.lvl[j,consumers], 0.50))],
                          size = rand_select[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[(troph.lvl[j,consumers] > quantile(troph.lvl[j,consumers], 0.50)) &
                                      (troph.lvl[j,consumers] < quantile(troph.lvl[j,consumers], 0.75))],
                          size = rand_select[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[(troph.lvl[j,consumers] > quantile(troph.lvl[j,consumers], 0.75))],
                          size = rand_select[k]/4))
    
    
    extinctions.rand_matTL[j,k] <- sum(extinctions[j, slct_rand]) #each 
    shannon.rand_matTL[j,k] <- diversity((abundances[j, slct_rand]))
  }
  
}


extinctions.rand_mat[1,1:9]

###
####4 DATA table with all response variables####
###
#test


#create df out with first connectance as first column, extinctions in the whole food web as second column
data <- cbind(con, as.vector(extinctions_mat)) %>%
  as.data.frame()
colnames(data) <- c("con", "ext_all")
head(data)

#add extinctions in sub-foodwebs
data <- cbind(data, 
              as.vector(n.bas_mat),
              as.vector(extinctions.big_mat),
              as.vector(extinctions.small_mat),
              as.vector(extinctions.BAS_mat))
colnames(data)[3:ncol(data)] <- c("n_basal" ,"ext_big", "ext_small", "ext_BAS")


#add shannon indices
data <- cbind(data,
              as.vector(shannon_mat),
              as.vector(shannon.big_mat),
              as.vector(shannon.small_mat),
              as.vector(shannon.BAS_mat))
colnames(data)[(ncol(data)-3):ncol(data)] <- c("shan_all", "shan_big", "shan_small", "shan_BAS")
head(data)


#add randomly drawn shannon samples
for (i in 1:length(rand_select)) {
  new <- shannon.rand_mat[,i] %>% as.vector()
  data[,ncol(data)+ 1 ] <- new
  colnames(data)[ncol(data)] <- c(paste("rand_shan", rand_select[i], sep = "_"))
}

#add randomly drawn extinction samples
for (i in 1:length(rand_select)) {
  new <- extinctions.rand_mat[,i] %>% as.vector()
  data[,ncol(data)+ 1 ] <- new
  colnames(data)[ncol(data)] <- c(paste("rand_ext", rand_select[i], sep = "_"))
}

#add randomly drawn TL considering shannon samples
for (i in 1:length(rand_select)) {
  new <- shannon.rand_matTL[,i] %>% as.vector()
  data[,ncol(data)+ 1 ] <- new
  colnames(data)[ncol(data)] <- c(paste("rand_shanTL", rand_select[i], sep = "_"))
}


#add randomly drawn TL considering extinction samples
for (i in 1:length(rand_select)) {
  new <- extinctions.rand_matTL[,i] %>% as.vector()
  data[,ncol(data)+ 1 ] <- new
  colnames(data)[ncol(data)] <- c(paste("rand_extTL", rand_select[i], sep = "_"))
}


####5 - quick check on correlation####  
cor(data$ext_all, data$rand_extTL_8) 
cor(data$ext_all, data$rand_ext_8)

cor(data$ext_all, data$ext_big)
cor(data$ext_all, data$ext_small)

cor(data$shan_all, data$rand_shan_32)
cor(data$shan_all, data$rand_shanTL_32)

cor(data$ext_all, data$ext_big)
cor(data$ext_all, data$ext_small)

####5 save data as .csv####
write.csv(data, "./data/20220505_96spec_15cons_v01.csv")





####Crap to go####
NULL
#actually,lets create a second set of points in the plots above


?TroLev()

#Neo Martinez

###
#model doesnt converge: change parameters until it convergence


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


