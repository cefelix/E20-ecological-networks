library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


###
####1 load raw RDS files, initialize parameters####
###

#1.1.1 read RDS list
d.list <- readRDS("./raw/output_64_2000randoms.rds")

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

#2.1.1 calculate extinction matrices
extinctions_mat <- apply(extinctions, MARGIN = c(1,2), FUN = sum) %>% #sum of extinctions for each connectance/replicate combination
  as.matrix()
extinctions.big_mat <- apply(extinctions[,,slct_bi], MARGIN = c(1,2), FUN = sum) %>% 
  as.matrix() 

extinctions.small_mat <- array(NA, dim = c(reps, length(con) ))
extinctions.BAS_mat <- array(NA, dim = c(reps, length(con) ))

for(i in 1:reps){
  for(j in 1:length(con)) {
    bas <- n_tot - sum(troph.lvl[i,j,] > 1) 
      #counts basal species in each food web
    extinctions.small_mat[i,j] <- sum(extinctions[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #selects extinctions in the smallest consumer species
    extinctions.BAS_mat[i,j] <- sum(extinctions[i,j, c(1:bas)])
      #selects extinctions in the basal species
  }
}


#2.2 Shannon indices
shannon_mat <-  apply(abundances, MARGIN = c(1,2), FUN = diversity) %>%
  as.matrix() #whole food web
shannon.big_mat <- apply(abundances[,,slct_bi], MARGIN = c(1,2), FUN = diversity) %>% 
  as.matrix() 

shannon.small_mat <- array(NA, dim = c(reps, length(con) ))
shannon.BAS_mat <- array(NA, dim = c(reps, length(con) ))

for(i in 1:reps){
  for(j in 1:length(con)) {
    bas <- n_tot - sum(troph.lvl[i,j,] > 1) 
    #counts basal species in each food web
    shannon.small_mat[i,j] <- diversity(abundances[i,j, c((bas+1):(bas+perc_sub*n_tot))])
      #calculates shannon index in the smallest consumer species
    shannon.BAS_mat[i,j]   <- diversity(abundances[i,j, c(1:bas)])
      #calculates shannon index in the basal species
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
              as.vector(extinctions.big_mat),
              as.vector(extinctions.small_mat),
              as.vector(extinctions.BAS_mat))
colnames(data)[3:ncol(data)] <- c("ext_big", "ext_small", "ext_BAS")


#add shannon indices
data <- cbind(data,
              as.vector(shannon_mat),
              as.vector(shannon.big_mat),
              as.vector(shannon.small_mat),
              as.vector(shannon.BAS_mat))
colnames(data)[(ncol(data)-3):ncol(data)] <- c("shan_all", "shan_big", "shan_small", "shan_BAS")
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
write.csv(data, "./data/64spec_2000randoms.csv")


####OLD####

#2.2 transform extinction output data to .csv

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

extinctions <- biomasses
extinctions[extinctions <= ext_thresh] <- 1
extinctions[extinctions > ext_thresh] <- 0
summary(extinctions) #not working yet



####4 - calculate SHANNON indices ####
####
####4.1 create subset foodwebs' abundance arrays####

abundances.BAS <- abundances[1:reps, 1:length(con), 1:n_bas]
abundances.small <- abundances[1:reps, 1:length(con), (n_bas+1):(n_bas+n_sub)]
abundances.big <- abundances[1:reps, 1:length(con), (n_tot-n_sub+1):n_tot]

####4.2 apply vegan::diversity over each abundance array, store output in .csv####
Sindex_all <- apply(abundances, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_BAS <- apply(abundances.BAS, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_big <- apply(abundances.small, MARGIN = c(1,2), FUN = diversity) %>%
  as.data.frame()
Sindex_small <- apply(abundances.big, MARGIN = c(1,2), FUN = diversity) %>%
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
  abuns_small = abundances[, con, (n_bas+1):(n_bas+n_sub)]
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


