library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


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


