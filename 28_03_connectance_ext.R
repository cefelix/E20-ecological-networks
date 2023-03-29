####0.0 libraries####

library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


####3.1 CONNECTANCE ON EXTINCTIONS/ABUNDANCES using a niche-model####
###
###
####3.3 set up parameters and food web####
n_bas <- 8                #no. of basal species
n_tot <- 128               #no. of all species

perc_sub <- 0.25          #share of the species present in sub-food webs (float between 0 and 1)
n_sub <- perc_sub*n_tot   #no. of species in a sub-foodweb

con <- seq(0.05, 0.45, by = 0.025)      #connectance of the food webs 
                                        #check out https://www.pnas.org/doi/epdf/10.1073/pnas.192407699
                                        #mininmal connectance 0.026, maximal 0.315 (from 17 empirical food webs, fig.1)

reps<- 1                             #no of replicates per food web
times <- seq(1, 1e10, by = 1e8)         #time for integration of dynamics
biom <- runif(n_tot, 1, 4)              #initial biomasses
BM <- runif(n_tot, 2, 3) %>%            #Body masses 
  sort() 
BM <- (10^BM)            


####3.4 initialize output objects, loop counts####
###
###

#output matrices: extinctions
extinctions_mat <- NULL         #extinctions in whole food web
extinctions.BAS_mat = NULL      #extinctions in the basal species
extinctions.small_mat = NULL    #extinctions in the 25% smallest species
extinctions.big_mat = NULL      #extinctions in the 25% biggest species

#output arrays: abundances (or: densities)

abundance_array <- array(dim = c(reps,length(con), n_tot))
abundance_array <- provideDimnames(abundance_array, sep = "_", base = list("rep", "con", "spec"))

#output matrices: shannon index
shannon_mat = NULL 

#loop counters (i set them up against the alphabet because that's cooler)
j=0
i=0


####3.5 compute 1e10 food webs####
for (j in 1:reps) {
  BM <- runif(n_tot, 2, 3) %>% #body masses of the species
    sort()
  BM <- (10^BM)
  extinctions <- rep(NA, length(con))
  
  extinctions.BAS <- rep(NA, length(con))
  extinctions.l <- rep(NA, length(con))
  extinctions.h <- rep(NA, length(con))
  
  shannon <- rep(NA, length(con))
  shannon_BAS <- rep(NA, length(con))
  shannon_small <- rep(NA, length(con))
  shannon_big <- rep(NA, length(con))
  
  i=0
  
  for (i in 1:length(con)){
    fw <- create_niche_model(S = n_tot, C = con[i])
    
    
    model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Unscaled()
    model$ext <- 0.1**6
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    
    #extinctions in entire food web
    exts = sol[nrow(sol), -1]
    exts = sum(exts <= model$ext)
    extinctions[i] = exts
    
    #extinctions in the basal species
    exts_BAS = sol[,-1]
    exts_BAS = exts_BAS[nrow(sol), 1:n_bas]
    exts_BAS = sum(exts_BAS <= model$ext)
    
    #smallest 16 species (without basals)
    exts_low = sol[,-1]
    slct = c((n_bas+1):(n_bas+n_sub))
    exts_low =exts_low[nrow(sol), slct]
    exts_low = sum(exts_low <= model$ext)
    
    #biggest 16 species
    exts_high = sol[,-1]
    slct = c((n_tot-n_sub+1):n_tot)
    exts_high = exts_high[nrow(sol), slct]
    exts_high = sum(exts_high <=model$ext)
    
    #storing the subset food webs in a vector
    extinctions.BAS[i] = exts_BAS
    extinctions.l[i] = exts_low
    extinctions.h[i] = exts_high
    
    #storing abundances in the abundance array
    abuns = sol[nrow(sol),-1] / BM
    abundance_array[j,i,] = abuns
    
    
    #storing shannon indices in the shannon vector 
    shannon[i] = diversity(abuns)
    
    
  }
  extinctions_mat = rbind(extinctions_mat, extinctions)
  
  extinctions.BAS_mat = rbind(extinctions.BAS_mat, extinctions.BAS)
  extinctions.small_mat = rbind(extinctions.small_mat, extinctions.l)
  extinctions.big_mat = rbind(extinctions.big_mat, extinctions.h)
  
  shannon_mat = rbind(shannon_mat, shannon)
  
  
  
  print(j)
  
}



####3.6 transform extinction output data to .csv####

extinctions_mat <- extinctions_mat %>% as.data.frame()
extinctions.big_mat <- extinctions.big_mat %>% as.data.frame()
extinctions.small_mat <- extinctions.small_mat %>% as.data.frame()
extinctions.BAS_mat <- extinctions.BAS_mat %>% as.data.frame()

colnames(extinctions_mat) <- paste0("Connectance_",con)
colnames(extinctions.big_mat) <- paste0("Connectance_",con)
colnames(extinctions.small_mat) <- paste0("Connectance_",con)
colnames(extinctions.BAS_mat) <- paste0("Connectance_",con)

getwd()
write.csv(extinctions_mat, "extinctions128all.csv")
write.csv(extinctions.big_mat, "extinctions128big.csv")
write.csv(extinctions.small_mat, "extinctions128small.csv")
write.csv(extinctions.BAS_mat, "extinctions128BAS.csv")

####3.7 transform abundance and shannon output data to .csv####
shannon_mat <- shannon_mat %>% as.data.frame()
colnames(shannon_mat) <- paste0("Connectance_", con)

write.csv(shannon_mat, "shannon_mat128.csv")

####4 calculate shannon indices for  abundance array####
####INSERT####

#"rep", "con", "spec"
for (k in 1:length(con)) {
  abuns_small = abundance_array[, con, (n_bas+1):(n_bas+n_sub)]
  apply(abuns_small, MARGIN = 1,diversity )
}

####5  plotting EXTINCTION OVER CONNECTANCE (entire food web)####


####5.1 load stored output#### 
#extinctions_mat <- read.csv("./extinctions100rep002.csv")[,-1]
#extinctions.big_mat <-read.csv("./extinctions100big.csv")
#extinctions.small_mat <- read.csv("./extinctions.small_mat")
#extinctions.BAS_mat <- read.csv("./extinctions100BAS.csv")

#create a df with connectance, extinctions and variance of extinctions
data.ext_con <- con     #set up df, first col will be connectance
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, mean))) #add mean no. of extinctions of n=rep iterations
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, var)))  #add variance of extinctions of n=rep iterations
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.big_mat, MARGIN = 2, mean))) #add mean no. of extinctions in biggest species
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.small_mat, MARGIN = 2, mean))) #add mean no. of extinctions in smallest species,
    #excluding Basal species
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.BAS_mat, MARGIN = 2, mean))) #add mean no. of extinctions in smallest species
    #including BASAL species
data.ext_con <- data.ext_con %>% #transpose and save as dataframe
  t() %>%
  as.data.frame()
  
#rename rows and columns:
rownames(data.ext_con) <- NULL
colnames(data.ext_con) <- c("con", "ext.m", "ext.v", "ext.big", "ext.small", "ext.BAS")

#calculate extinction rates
data.ext_con$rate <- data.ext_con$ext.m/n_tot
data.ext_con$ext.big <- data.ext_con$ext.big/n_sub
data.ext_con$ext.small <- data.ext_con$ext.small/n_sub
data.ext_con$ext.BAS <- data.ext_con$ext.BAS/n_sub
  
data.ext_con %>% str()
  
#plot the mean extinction numbers
ggplot(data.ext_con, aes(x = con, y=ext.m))+
  geom_point()+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.m)+1)) #becomes asymptotic to x-axis parallel with connectance -> 0.5

#plot the variance of the extinctions
ggplot(data.ext_con, aes(x = con, y=ext.v))+
  geom_point()+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.v)+1)) #interesting pattern, lets have a look at it with 10000 rep's

#plot the rate of extinctions
ggplot(data.ext_con, aes(x = con, y=rate) )+
  geom_point()+
  geom_point(aes(x=con, y=ext.BAS, color = "basal species"))+
  geom_point(aes(x=con, y=ext.small, color = "smallest consumers"))+
  geom_point(aes(x=con, y=ext.big, color = "biggest consumers"))+
  #geom_point(aes(x=con, y=ext.BAS, color = "red"))+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.small)+0.1)) #becomes asymptotic to x-axis parallel with connectance -> 0.5
  


####6 correlations between whole food web and sub food webs####  






####666####
NULL
#actually,lets create a second set of points in the plots above


?TroLev()

#Neo Martinez


