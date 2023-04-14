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

con <- seq(0.05, 0.45, by = 0.025)      #connectance of the food webs 
                                        #check out https://www.pnas.org/doi/epdf/10.1073/pnas.192407699
                                        #mininmal connectance 0.026, maximal 0.315 (from 17 empirical food webs, fig.1)

reps<- 1                             #no of replicates per food web
times <- seq(1, 1e12, by = 1e9)         #time for integration of dynamics
biom <- runif(n_tot, 1, 4)              #initial biomasses
BM <- runif(n_tot, 1, 12) %>%           #Body masses 
  sort()                                #realistic spanning range in soil food webs (Potapov 2021):
BM <- (10^BM)                           #https://pubmed.ncbi.nlm.nih.gov/34086977/
ext_thresh <- 0.1**6                     #threshold below which species is considered, extinct


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
saveRDS(abundance_array, file = "./raw/abundance.rds")
saveRDS(biomass_array, file = "./raw/biomass.rds")
saveRDS(extinction_array, file = "./raw/extinction.rds")
saveRDS(troph.lvl_array, file = "./raw/trophic_lvl.rds")


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


###
####3 calculate EXTINCTION matrices####
###

#3.1.0 initialize selection variables
slct_bi <- c((n_tot-n_sub+1):n_tot)       #selects biggest consumer species (n_sub=32) 
slct_sm <- c((n_bas+1):(n_sub+n_bas))     #selects smallest consumer species (n_sub=32) 
slct_BAS <- c(1:n_bas)                    #selects basal species

#3.1.1 calculate extinction matrices
extinctions_mat <- apply(extinction_array, MARGIN = c(1,2), FUN = sum) %>% #sum of extinctions for each connectance/replicate combination
  as.data.frame()  

extinctions.big_mat <- apply(extinction_array[,,slct_bi], MARGIN = c(1,2), FUN = sum) %>% 
  as.data.frame() 
extinctions.small_mat <- apply(extinction_array[,,slct_sm], MARGIN = c(1,2), FUN = sum) %>% 
  as.data.frame() 
extinctions.BAS_mat <- apply(extinction_array[,,slct_BAS], MARGIN = c(1,2), FUN = sum) %>% 
  as.data.frame() 

  #apply(extinction_array, MARGIN = 3, FUN = sum) #(sum of extinctions for each species)


#3.1.2 transform extinction output data to .csv

colnames(extinctions_mat) <- paste0("Connectance_",con)
colnames(extinctions.big_mat) <- paste0("Connectance_",con)
colnames(extinctions.small_mat) <- paste0("Connectance_",con)
colnames(extinctions.BAS_mat) <- paste0("Connectance_",con)

write.csv(extinctions_mat, "./data/extinctions128all.csv")
write.csv(extinctions.big_mat, "./data/extinctions128big.csv")
write.csv(extinctions.small_mat, "./data/extinctions128small.csv")
write.csv(extinctions.BAS_mat, "./data/extinctions128BAS.csv")


####3.7 transform abundance and shannon output data to .csv####
shannon_mat <- shannon_mat %>% as.data.frame()
colnames(shannon_mat) <- paste0("Connectance_", con)

write.csv(shannon_mat, "./real/shannon_mat128.csv")


####3.8 transform extinction output to relevant subset data frames

extinction_array <- biomass_array
extinction_array[extinction_array <= ext_thresh] <- 1
extinction_array[extinction_array > ext_thresh] <- 0
summary(extinction_array) #not working yet



####4 - calculate SHANNON indices ####
####
####4.1 create subset foodwebs' abundance arrays####

abundance_array.BAS <- abundance_array[1:reps, 1:length(con), 1:n_bas]
abundance_array.small <- abundance_array[1:reps, 1:length(con), (n_bas+1):(n_bas+n_sub)]
abundance_array.big <- abundance_array[1:reps, 1:length(con), (n_tot-n_sub+1):n_tot]

####4.2 apply vegan::diversity over each abundance array, store output in .csv####
Sindex_all <- apply(abundance_array, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_BAS <- apply(abundance_array.BAS, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_big <- apply(abundance_array.small, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_small <- apply(abundance_array.big, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()

#write as csv:
write.csv(Sindex_all, "./real/Sindex_all.csv")
write.csv(Sindex_BAS, "./real/Sindex_BAS.csv")
write.csv(Sindex_big, "./real/Sindex_big.csv")
write.csv(Sindex_small, "./real/Sindex_small.csv")


####5  plotting EXTINCTION OVER CONNECTANCE (entire food web)####


####5.1 load stored output#### 

extinctions_mat <- read.csv("./extinctions128all.csv")[,-1]
#extinctions.big_mat <-read.csv("./extinctions128big.csv")[,-1]
#extinctions.small_mat <- read.csv("./extinctions128small.csv")[,-1]
#extinctions.BAS_mat <- read.csv("./extinctions128BAS.csv")[,-1]

#shannon_mat <- read.csv("./shannon_mat128.csv")[,-1]

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

####6 plotting shannon index over connectance####
data.shan <- con

data.shan <- rbind(data.shan, t(apply(Sindex_all, MARGIN = 2, mean)))
data.shan <- rbind(data.shan, t(apply(Sindex_BAS, MARGIN = 2, mean)))
data.shan <- rbind(data.shan, t(apply(Sindex_big, MARGIN = 2, mean)))
data.shan <- rbind(data.shan, t(apply(Sindex_small, MARGIN = 2, mean)))
data.shan <- data.shan %>% #transpose and save as dataframe
  t() %>%
  as.data.frame()

colnames(data.shan) <- c("con", "all", "basal", "big25", "small25")

#plot it:

ggplot(data.shan, aes(x = con, y=all) )+
  geom_point()+
  geom_point(aes(x=con, y=basal, color = "basal species"))+
  geom_point(aes(x=con, y=small25, color = "smallest consumers"))+
  geom_point(aes(x=con, y=big25, color = "biggest consumers"))+
  #geom_point(aes(x=con, y=ext.BAS, color = "red"))+
  scale_y_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, log(n_tot)+1))+
  ylab("Shannon Index")+
  xlab("Connectance")
  #geom_abline(intercept = log(128), slope = 0, aes(color="black"))+
  #geom_abline(intercept = log(32), slope = 0, aes(color=F8766D))
  


####7 correlations between whole food web and sub food webs####  


cor(data.shan$all, data.shan$basal)
cor(data.shan$all, data.shan$big25)
cor(data.shan$all, data.shan$small25)

cor(data.ext_con$rate, data.ext_con$ext.BAS) #Na, because  no extinctions in the basal species
cor(data.ext_con$rate, data.ext_con$ext.big)
cor(data.ext_con$rate, data.ext_con$ext.small)




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


