library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


###
####1 load raw RDS files, initialize parameters####
###

#1.1.1 read RDS list
d.list <- readRDS("./raw/output_96_2000randoms.rds")

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
#slct_sm <- c((n_bas+1):(n_sub+n_bas))    #OUTDATED: selects smallest consumer species (n_sub=32) 
#slct_BAS <- c(1:n_bas)                   #OUTDATED: selects basal species



###
####2 response variable matrices####
### 

#2.1 initializing output arrays
n.bas_mat <- array(NA, dim = c(reps, length(con) ))


extinctions_mat <- apply(extinctions, MARGIN = c(1,2), FUN = sum) %>% #sum of extinctions for each connectance/replicate combination
  as.matrix()
extinctions.small_mat <- array(NA, dim = c(reps, length(con) ))
extinctions.BAS_mat <- array(NA, dim = c(reps, length(con) ))
extinctions.big_mat <- array(NA, dim = c(reps, length(con) ))
extinctions.rand_mat <- array(NA, dim = c(reps, length(con) ))

for(i in 1:reps){
  for(j in 1:length(con)) {
    bas <- n_tot - sum(troph.lvl[i,j,] > 1) 
      #counts basal species in each food web
    extinctions.small_mat[i,j] <- sum(extinctions[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #selects extinctions in the smallest consumer species
    extinctions.BAS_mat[i,j] <- sum(extinctions[i,j, c(1:bas)])
      #selects extinctions in the basal species
    extinctions.big_mat[i,j] <- sum(extinctions[i,j, slct_bi])
      #extinctions in the biggest consumer species
    
    consumers <- c((bas+1):n_tot)                   # a vector of all consumer species
    rand_sub <- sample(consumers, size = n_sub)     # a random selection of n_sub consumer species
    extinctions.rand_mat[i,j] <- sum(extinctions[i,j, rand_sub]) #
  }
}


shannon_mat <-  apply(abundances, MARGIN = c(1,2), FUN = diversity) %>%
  as.matrix() #whole food web
shannon.small_mat <- array(NA, dim = c(reps, length(con) ))
shannon.BAS_mat <- array(NA, dim = c(reps, length(con) ))
shannon.big_mat <- array(NA, dim = c(reps, length(con) ))
shannon.rand_mat <- array(NA, dim = c(reps, length(con) )) 


#2.2 calculating extinctions and shannon indices
for(i in 1:reps){
  for(j in 1:length(con)) {
    
    bas <- n_tot - sum(troph.lvl[i,j,] > 1)         # a vector of all basal species
    consumers <- c((bas+1):n_tot)                   # a vector of all consumer species
    slct_rand <- sample(consumers, size = n_sub)    # a random selection of consumer species
    n.bas_mat[i,j] <- bas                           # the number of basal species in the system
    
   
    extinctions.small_mat[i,j] <- sum(extinctions[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      # selects extinctions in the smallest consumer species
    extinctions.BAS_mat[i,j] <- sum(extinctions[i,j, c(1:bas)])
      # selects extinctions in the basal species
    extinctions.big_mat[i,j] <- sum(extinctions[i,j, slct_bi])
      # extinctions in the biggest consumer species
    extinctions.rand_mat[i,j] <- sum(extinctions[i,j, slct_rand]) 
      # extinctions in a random selection of consumer species
    
    
    shannon.small_mat[i,j] <- diversity(abundances[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      # calculates shannon index in the smallest consumer species
    shannon.BAS_mat[i,j]   <- diversity(abundances[i,j, c(1:bas)])

      # calculates shannon index in the basal species
    shannon.big_mat[i,j] <- diversity(abundances[i,j, slct_bi])
      # shannon index in the biggest consumer species
    shannon.rand_mat[i,j] <- diversity((abundances[i,j, slct_rand]))

  }
}




#2.3 Trophic levels
  troph.max_mat <- apply(troph.lvl, MARGIN = c(1,2), FUN = max) %>%
    as.matrix()
  troph.min_mat <- apply(troph.lvl, MARGIN = c(1,2), FUN = min) %>%
    as.matrix()
  troph.var_mat <- apply(troph.lvl, MARGIN = c(1,2), FUN = var) %>%
    as.matrix()
  
  troph.big.max_mat <- apply(troph.lvl[,,slct_bi], MARGIN = c(1,2), FUN = max) %>% 
    as.matrix() 
  troph.big.min_mat <- apply(troph.lvl[,,slct_bi], MARGIN = c(1,2), FUN = min) %>% 
    as.matrix() 
  troph.big.var_mat <- apply(troph.lvl[,,slct_bi], MARGIN = c(1,2), FUN = var) %>% 
    as.matrix()
  
  troph.small.max_mat <- array(NA, dim = c(reps, length(con)))
  troph.small.min_mat <- array(NA, dim = c(reps, length(con)))
  troph.small.var_mat <- array(NA, dim = c(reps, length(con)))
  
for(i in 1:reps){
  for(j in 1:length(con)) {
    bas <- n_tot - sum(troph.lvl[i,j,] > 1) 
      #counts basal species in each food web
    troph.small.max_mat[i,j] <- max(troph.lvl[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #selects max trophic levels in the smallest consumer species
    troph.small.min_mat[i,j] <- min(troph.lvl[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #min()
    troph.small.var_mat[i,j] <- var(troph.lvl[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #var()
  }
}  


###
####4 DATA table with all response variables####
###
#test

#create df out with first connectance as first column, extinctions in the whole food web as second column
data <- cbind(sort(rep(con, reps)), as.vector(extinctions_mat)) %>%
  as.data.frame()
colnames(data) <- c("con", "ext_all")
  head(data)
  
#add extinctions in sub-foodwebs
data <- cbind(data, 
              as.vector(n.bas_mat),
              as.vector(extinctions.big_mat),
              as.vector(extinctions.small_mat),
              as.vector(extinctions.BAS_mat),
              as.vector(extinctions.rand_mat))
colnames(data)[3:ncol(data)] <- c("n_basal" ,"ext_big", "ext_small", "ext_BAS", "ext_rand")


#add shannon indices
data <- cbind(data,
              as.vector(shannon_mat),
              as.vector(shannon.big_mat),
              as.vector(shannon.small_mat),
              as.vector(shannon.BAS_mat),
              as.vector(shannon.rand_mat))
colnames(data)[(ncol(data)-4):ncol(data)] <- c("shan_all", "shan_big", "shan_small", "shan_BAS", "shan_rand")
  head(data)

#add trophic level information
data <- cbind(data,
               as.vector(troph.max_mat),
               as.vector(troph.min_mat),
               as.vector(troph.var_mat),
               
               as.vector(troph.big.max_mat),
               as.vector(troph.big.min_mat),
               as.vector(troph.big.var_mat),
               
               as.vector(troph.small.max_mat),
               as.vector(troph.small.min_mat),
               as.vector(troph.small.var_mat)
)

#add colnames
colnames(data)[(ncol(data)-8):ncol(data)] <- c(
  "maxTL", "minTL", "varTL",
  "maxTL_B", "minTL_B", "varTL_B",
  "maxTL_s", "minTL_s", "varTL_s"
)



####5 save data as .csv####
write.csv(data, "./data/96spec_2000randoms.csv")


####666####
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


