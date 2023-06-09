library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


###
####1 load raw RDS files, initialize parameters####
###

#1.1.1 read RDS list
list.files("./raw")
d.list <- readRDS("./raw/20220506_96spec_2000c_04to12bas_v01.rds")


  abundances <- d.list$abundances
  biomasses <- d.list$biomasses
  extinctions <- d.list$extinctions
  troph.lvl <- d.list$troph_lvl
  
  food_START <- d.list$prey_on_START
  food_END <- d.list$prey_on_END
  consumed_bySTART <- d.list$consumed_by_START
  consumed_byEND <- d.list$consumed_by_END


#1.2 initialize selection variables


#first initialize model parameters: go to  ODE-script
slct_bi <- c((n_tot-n_sub+1):n_tot)       #selects biggest consumer species (n_sub=32) 


###
####2 response variable matrices####
### 

#2.1 initializing output arrays
extinctions_vec <- apply(extinctions, MARGIN = c(1), FUN = sum) %>% #sum of extinctions for each connectance/replicate combination
  as.vector()
  extinctions.small_vec   <- rep(NA, length(con))
  extinctions.BAS_vec     <- rep(NA, length(con))
  extinctions.big_vec     <- rep(NA, length(con))

shannon_vec <-  apply(abundances, MARGIN = c(1), FUN = diversity) %>%
  as.vector() #whole food web
  shannon.small_vec   <- rep(NA, length(con))
  shannon.BAS_vec     <- rep(NA, length(con))
  shannon.big_vec     <- rep(NA, length(con))
  
n.bas_vec <- rep(NA, length(con))  


#2.2 calculating extinctions and shannon indices
  for(j in 1:length(con)) {
    
    bas <- n_tot - sum(troph.lvl[j,] > 1)             # the amount of basal species
    consumers <- c((bas+1):n_tot)                     # a vector of all consumer species
    n.bas_vec[j] <- bas                               # the number of basal species in the system
    
    slct_sm <- consumers[1:n_sub]                     # the smallest consumer species
    slct_bi <- c((n_tot-n_sub+1):n_tot)               # the biggest consumer species
    
    
  #extinctions in 24 species
    extinctions.small_vec[j] <- sum(extinctions[j, slct_sm])
      # selects extinctions in the smallest consumer species
    extinctions.BAS_vec[j] <- sum(extinctions[j, c(1:bas)])
      # selects extinctions in the basal species
    extinctions.big_vec[j] <- sum(extinctions[j, slct_bi])
      # extinctions in the biggest consumer species

  #diversity in 24 species
    shannon.small_vec[j] <- diversity(abundances[j, slct_sm])
      # calculates shannon index in the smallest consumer species
    shannon.BAS_vec[j]   <- diversity(abundances[j, c(1:bas)])
      # calculates shannon index in the basal species
    shannon.big_vec[j] <- diversity(abundances[j, slct_bi])
      # shannon index in the biggest consumer species
    
    
    
    
  }


#output for random selection

sample_size <- seq(from = 36, to = 4, by=-4) #the sample sizes for the random selections
shannon.rand_mat        <- array(NA, dim = c( length(con), length(sample_size) )) #randomly selected 
extinctions.rand_mat    <- array(NA, dim = c( length(con), length(sample_size) ))

shannon.rand_matTL      <- array(NA, dim = c( length(con), length(sample_size) )) #randomly selected, but from different trophic levels
extinctions.rand_matTL  <- array(NA, dim = c( length(con), length(sample_size) ))

for (j  in 1:length(con)) {
  bas <- n_tot - sum(troph.lvl[j,] > 1)       # the number of basal species
  consumers <- c((bas+1):n_tot)               # a vector of all consumer species
  n.bas_vec[j] <- bas
  
  #random selection of species:
  for (k in 1:length(sample_size)) {
    slct_rand <- sample(consumers, size = sample_size[k]) %>%
      sort() # a random selection of consumer species
    
    extinctions.rand_mat[j,k] <- sum(extinctions[j, slct_rand])     
      # sum of extinctions at con=j for k selections of random species
    shannon.rand_mat[j,k] <- diversity((abundances[j, slct_rand]))  
      #shannon index at con=j
  }
  
  #random selection of species, but from different trophic levels:
  TL_c <- troph.lvl[j,consumers] 
  
  for (k in 1:length(sample_size)) {
    
    slct_rand <-  sample(consumers[TL_c <= quantile(TL_c, 0.25)],
                          size = sample_size[k]/4)
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.25) &
                                    TL_c <= quantile(TL_c, 0.5)],
                          size = sample_size[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.5) &
                                    TL_c <= quantile(TL_c, 0.75)],
                          size = sample_size[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.75)],
                          size = sample_size[k]/4))
    
    extinctions.rand_matTL[j,k] <- sum(extinctions[j, slct_rand]) #each 
    shannon.rand_matTL[j,k] <- diversity((abundances[j, slct_rand]))
  }
  print(j)
}

####3.1 compare exticntion rates/ realized shannon indices for each food web####

rate_comp.bi <- rep(NA, length(con))
rate_comp.sm <- rep(NA, length(con))

shan_comp.bi <- rep(NA, length(con))
shan_comp.sm <- rep(NA, length(con))

for (i in 1:length(con)) {
  n_bas <- sum(troph.lvl[i]== 1)
  exts_all <- extinctions_vec[i]/(96-n_bas) #true ext rate in the total sample
  exts_big <- extinctions.big_vec[i]/24     #ext rate bigs
  exts_sml <- extinctions.small_vec[i]/24   #ext rate smalls
  
  rate_comp.bi[i] <- exts_big / exts_all
  rate_comp.sm[i] <- exts_sml / exts_all
  
  shan_all      <- shannon_vec[i]
  shan_big  <- shannon.big_vec[i]
  shan_sml  <- shannon.small_vec[i]
  
  shan_comp.bi[i] <- shan_big/shan_all
  shan_comp.sm[i] <- shan_sml/shan_all
  
}



plot(con, rate_comp.bi) #sampling from big species overestimates extinction rates
plot(con, rate_comp.sm) #sampling from small species underestimates extinction rates

plot(con, shan_comp.bi, ylim = c(0, 2))
plot(con, shan_comp.sm, ylim = c(0, 2))

####3.2 for varying sub-sample size####

ext_comps <- array(data=NA, 
                   dim = c(length(con), length(sample_size), 2) ) #3rd dimension indicates whether trophic levels were considered

shan_comps <- array(data=NA, 
                    dim = c(length(con), length(sample_size), 2) )

for (i in 1:length(con)) {
  n_bas <- sum(troph.lvl[i]== 1)
  exts.ALL <- extinctions_vec[i]/(96-n_bas)
  shans_all <- shannon_vec[i]/log(96) #Shannon equitability index 
  
  for (j in 1:length(sample_size)){
    exts    <- extinctions.rand_mat[i, j]/sample_size[j]
    extsTL  <- extinctions.rand_matTL[i, j]/sample_size[j]
    
    ext_comps[i, j, 1] <- exts/exts.ALL
    ext_comps[i, j, 2] <- extsTL/exts.ALL
    
    ####ADD####
    
    shans   <- shannon.rand_mat[i,  j]/log(sample_size[j]-extinctions.rand_mat[i, j])
    shansTL <- shannon.rand_matTL[i,j]/log(sample_size[j]-extinctions.rand_matTL[i, j])   
      
    shan_comps[i, j, 1] <- shans/shans_all
    shan_comps[i, j, 2] <- shansTL/shans_all
  }
}


library(tidyr)
#saving output:
#extinction-rate ratios for random sampling
d.rand_exts   <- ext_comps[,,1] %>% as.data.frame()
d.rand_exts$con <- con
  colnames(d.rand_exts) <- c("36", "32", "28", "24", "20", 
                             "16", "12", "8", "4")
data <- gather(d.rand_exts, 
                      key= subsize, value=ext_ratio, 
                      "36":"4", factor_key=TRUE)

#extinction-rate ratios for TL sampling
d.rand_extsTL <- ext_comps[,,2] %>% as.data.frame()
colnames(d.rand_extsTL) <- c("36", "32", "28", "24", "20", 
                             "16", "12", "8", "4")
extsTL <- gather(d.rand_extsTL, 
                key= subsize, value=ext_ratioTL, 
                "36":"4", factor_key=TRUE)
data$ext_ratioTL <- extsTL$ext_ratioTL


#shannon indices from random sampling
shan.df1 <- shan_comps[,,1] %>% as.data.frame()
  colnames(shan.df1) <- c("36", "32", "28", "24", "20", 
                         "16", "12", "8", "4")
shan.df1 <- gather(shan.df1,
                   key=subsize, value= shannon_random,
                   "36":"4", factor_key = TRUE)  
data$shan_rand <- shan.df1$shannon_random


#shannon indices from TL-restricted sampling
shan.df2 <- shan_comps[,,2] %>% as.data.frame()
  colnames(shan.df2) <- c("36", "32", "28", "24", "20", 
                        "16", "12", "8", "4")
shan.df2 <- gather(shan.df2,
                   key=subsize, value= shannon_random,
                   "36":"4", factor_key = TRUE)   
data$shan_TL <- shan.df2$shannon_random

#check data
colnames(data)[1] <- "con"
str(data)
write.csv(data, "./data/20230523_SampleSizes_v01.csv")

###
####4 DATA table with all response variables####
###
#test


#create df out with first connectance as first column, extinctions in the whole food web as second column
data <- cbind(con, extinctions_vec) %>%
  as.data.frame()
colnames(data) <- c("con", "ext_all")
  head(data)
  
#add extinctions in sub-foodwebs
data <- cbind(data,
              n.bas_vec,
              extinctions.big_vec,
              extinctions.small_vec,
              extinctions.BAS_vec)
colnames(data)[3:ncol(data)] <- c("n_basal" ,"ext_big", "ext_small", "ext_BAS")


#add shannon indices
data <- cbind(data,
              shannon_vec,
              shannon.big_vec,
              shannon.small_vec,
              shannon.BAS_vec)
colnames(data)[(ncol(data)-3):ncol(data)] <- c("shan_all", "shan_big", "shan_small", "shan_BAS")
  head(data)

    
#add randomly drawn shannon samples
  for (i in 1:length(sample_size)) {
    new <- shannon.rand_mat[,i] %>% as.vector()
    data[,ncol(data)+ 1 ] <- new
    colnames(data)[ncol(data)] <- c(paste("rand_shan", sample_size[i], sep = "_"))
  }
  
#add randomly drawn extinction samples
  for (i in 1:length(sample_size)) {
    new <- extinctions.rand_mat[,i] %>% as.vector()
    data[,ncol(data)+ 1 ] <- new
    colnames(data)[ncol(data)] <- c(paste("rand_ext", sample_size[i], sep = "_"))
  }
  
#add randomly drawn TL considering shannon samples
  for (i in 1:length(sample_size)) {
    new <- shannon.rand_matTL[,i] %>% as.vector()
    data[,ncol(data)+ 1 ] <- new
    colnames(data)[ncol(data)] <- c(paste("rand_shanTL", sample_size[i], sep = "_"))
  }
  
  
#add randomly drawn TL considering extinction samples
  for (i in 1:length(sample_size)) {
    new <- extinctions.rand_matTL[,i] %>% as.vector()
    data[,ncol(data)+ 1 ] <- new
    colnames(data)[ncol(data)] <- c(paste("rand_extTL", sample_size[i], sep = "_"))
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
write.csv(data, "./data/20220506_96spec_2000c_04to12bas_v02.csv")


  
  
  
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


