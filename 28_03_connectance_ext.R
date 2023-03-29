####0.0 libraries####

library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


####3.1 CONNECTANCE ON EXTINCTIONS/ABUNDANCES using a niche-model####


####3.2 defining explanatory/response variables####  
#explanatory:
#connectance of the food web
#n_sub / n_tot = size of the sub-food web in relation to the total food web

#responses (each to be compared, between a food web of n species and a sub-foodweb of k<n species)
No.ext = NULL #number of extinctions
abundances =NULL #total biomass / Body mass

####3.3 set up parameters and food web####
n_bas <- 8
n_tot <- 64



#seeds.fw <- seq(1,length(con)*length(reps), by=1) #seeds for creating the food web, 
#LEAVE OUT FOR NOW

con <- seq(0.05, 0.45, by = 0.025) #connectance values of the food webs (EXPLANATORY)
reps<- 100 #no of replicates per food web
times <- seq(1, 1e8, by = 1e6) #times for integration
biom <- runif(n_tot, 1, 4) #initial biomasses - changing them doesn't affect the abundances in the long run
BM <- runif(n_tot, 2, 3) %>% #body masses of the species
  sort()
BM <- (10^BM)#does this affect the abundances?



#intialise variables for loop
i = 0
exts = c()

#initialize output vectors for whole food web:
extinctions <- rep(NA, length(con))

extinctions.lBAS <- rep(NA, length(con))
extinctions.l <- rep(NA, length(con))
extinctions.l <- rep(NA, length(con))

biomasses_out <- matrix(data=NA, nrow = n_tot, ncol = length(con)) 
  #a matrix containing the biomasses at the final time 
  #for each connectance value (columns) and species (rows)
abundances_out <- matrix(data = NA, nrow = n_tot, ncol = length(con))
extinctions_mat <- NULL

#for the sub-foodwebs
extinctions.lBAS_mat = NULL
extinctions.l_mat = NULL
extinctions.h_mat = NULL


#initialize steps to select a number for a subfoodweb:
#by = 0.25 #in relation to n_tot, should produce an integer when multiplying with n_tot
#n_sub = seq(by*n_tot, 48, by = by*n_tot) lets leave that for now and do it stupid but easy
n_sub = 16


#create a sub-foodweb, using the 25% lowest BM species:
#tresh = 0.25
#tresh = tresh*length(BM)

#sub.low <- BM[1:tresh]


#create a sub-foodweb, using the 25% upper BM species:
#tresh = 0.75
#tresh = tresh*length(BM)
#tresh = tresh+1 #these 3 rows are ugly, but otherwise it does not work -.-

#sub.high <- BM[tresh:length(BM)]
#NOW THIS HAS TO BE APPLIED ON SOL INSIDE THE LOOP!


#this line just creates a niche model with C=0.15
#fw <- create_niche_model(S = n_tot, C = 0.45)
i=0
j=0

for (j in 1:reps) {
  BM <- runif(n_tot, 2, 3) %>% #body masses of the species
    sort()
  BM <- (10^BM)
  extinctions <- rep(NA, length(con))
  
  extinctions.lBAS <- rep(NA, length(con))
  extinctions.l <- rep(NA, length(con))
  extinctions.h <- rep(NA, length(con))
  
  i=0
  
  for (i in 1:length(con)){
    fw <- create_niche_model(S = n_tot, C = con[i])
    
    
    model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Unscaled()
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    exts = sol[nrow(sol), -1]
    exts = sum(exts <= 0.1)
    extinctions[i] = exts
    print(i)
    #biomasses_out[,i] = sol[nrow(sol), -1]
    #abundances_out[,i] = sol[nrow(sol), -1]/BM
    #takes the biomasses per species and puts them into a matrix
    #shannon <- abundances_out %>%
      #round() %>%
      #diversity()
    
    ####
    #### CREATE A FOR-LOOP HERE, to compare a gradient of sub-foodwebs to whole food web
    ####
    
    ##look at smallest 16 species, including basal species
    exts_lowBAS = sol[,-1]
    exts_lowBAS = exts_lowBAS[nrow(sol), 1:n_sub]
    exts_lowBAS = sum(exts_lowBAS <= 0.1)
    
    #lowest 16 without basals
    exts_low = sol[,-1]
    slct = c((n_bas+1):(n_bas+n_sub))
    exts_low =exts_low[nrow(sol), slct]
    exts_low = sum(exts_low <= 0.1)
    
    
    print(i)
    
    #biggest 16
    exts_high = sol[,-1]
    slct = c((n_tot-n_sub+1):n_tot)
    exts_high = exts_high[nrow(sol), slct]
    exts_high = sum(exts_high <= 0.1)
    
    #storing the subset food webs in a vector
    extinctions.lBAS[i] = exts_lowBAS
    extinctions.l[i] = exts_low
    extinctions.h[i] = exts_high
    
    
    
  }
  extinctions_mat = rbind(extinctions_mat, extinctions)
  
  extinctions.lBAS_mat = rbind(extinctions.lBAS_mat, extinctions.lBAS)
  extinctions.l_mat = rbind(extinctions.l_mat, extinctions.l)
  extinctions.h_mat = rbind(extinctions.h_mat, extinctions.h)
  
  
  
  print(j)
  
}
extinctions_mat <- extinctions_mat %>% as.data.frame()
extinctions.h_mat <- extinctions.h_mat %>% as.data.frame()
extinctions.l_mat <- extinctions.l_mat %>% as.data.frame()
extinctions.lBAS_mat <- extinctions.lBAS_mat %>% as.data.frame()

colnames(extinctions_mat) <- paste0("Connectance_",con)
colnames(extinctions.h_mat) <- paste0("Connectance_",con)
colnames(extinctions.l_mat) <- paste0("Connectance_",con)
colnames(extinctions.lBAS_mat) <- paste0("Connectance_",con)

getwd()
write.csv(extinctions_mat, "extinctions100all.csv")
write.csv(extinctions.h_mat, "extinctions100big.csv")
write.csv(extinctions.l_mat, "extinctions100small.csv")
write.csv(extinctions.lBAS_mat, "extinctions100BAS.csv")


####4.1 plotting EXTINCTION OVER CONNECTENCE (entire food web)####
#extinctions_mat <- read.csv("./extinctions100rep002.csv")[,-1]

#create a df with connectance, extinctions and variance of extinctions
data.ext_con <- con     #set up df, first col will be connectance
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, mean))) #add mean no. of extinctions of n=rep iterations
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, var)))  #add variance of extinctions of n=rep iterations
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.h_mat, MARGIN = 2, mean))) #add mean no. of extinctions in biggest species
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.l_mat, MARGIN = 2, mean))) #add mean no. of extinctions in smallest species,
    #excluding Basal species
data.ext_con <- rbind(data.ext_con, t(apply(extinctions.lBAS_mat, MARGIN = 2, mean))) #add mean no. of extinctions in smallest species
    #including BASAL species
data.ext_con <- data.ext_con %>% #transpose and save as dataframe
  t() %>%
  as.data.frame()
  
#rename rows and columns:
rownames(data.ext_con) <- NULL
colnames(data.ext_con) <- c("con", "ext.m", "ext.v", "ext.big", "ext.small", "ext.smallBAS")

#calculate extinction rates
#data.ext_con$rate <- data.ext_con$ext.m/n_tot
data.ext_con$rate <- data.ext_con$ext.m/n_tot
data.ext_con$ext.big <- data.ext_con$ext.big/n_sub
data.ext_con$ext.small <- data.ext_con$ext.small/n_sub
data.ext_con$ext.smallBAS <- data.ext_con$ext.smallBAS/n_sub
  
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
  geom_point(aes(x=con, y=ext.smallBAS, color = "smallest including Basals"))+
  geom_point(aes(x=con, y=ext.small, color = "smallest consumers"))+
  geom_point(aes(x=con, y=ext.big, color = "biggest consumers"))+
  #geom_point(aes(x=con, y=ext.smallBAS, color = "red"))+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.small)+0.1)) #becomes asymptotic to x-axis parallel with connectance -> 0.5
  
  #this is partly nonsense, because i set an fixed extinction treshold, and it would have to be dependent on BM of species



####4.2 plotting the same stuff for a subset food web####
NULL
#actually,lets create a second set of points in the plots above







