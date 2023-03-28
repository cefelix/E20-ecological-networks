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
n_sub = seq(0,1, by=0.125)*(n_tot-n_bas) 


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

#initialize output vector:
extinctions <- rep(NA, length(con))
biomasses_out <- matrix(data=NA, nrow = n_tot, ncol = length(con)) 
  #a matrix containing the biomasses at the final time 
  #for each connectance value (columns) and species (rows)
abundances_out <- matrix(data = NA, nrow = n_tot, ncol = length(con))
extinctions_mat <- NULL
 


#this line just creates a niche model with C=0.15
#fw <- create_niche_model(S = n_tot, C = 0.45)
i=0
j=0

for (j in 1:reps) {
  BM <- runif(n_tot, 2, 3) %>% #body masses of the species
    sort()
  BM <- (10^BM)
  
  for (i in 1:length(con)){
    fw <- create_niche_model(S = n_tot, C = con[i])
    
    
    model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
      initialise_default_Unscaled()
    
    #solve ode's
    sol = lsoda_wrapper(times, biom, model)
    exts = sol[nrow(sol), -1]
    exts = sum(exts <= 0.1)
    extinctions[i] = exts
    #biomasses_out[,i] = sol[nrow(sol), -1]
    #abundances_out[,i] = sol[nrow(sol), -1]/BM
    #takes the biomasses per species and puts them into a matrix
    #shannon <- abundances_out %>%
      #round() %>%
      #diversity()
  }
  extinctions_mat = rbind(extinctions_mat, extinctions)
  print(j)
  
}
extinctions_mat <- extinctions_mat %>% as.data.frame()


colnames(extinctions_mat) <- paste0("Connectance_",con)


####4.1 plotting EXTINCTION OVER CONNECTENCE (entire food web)####

#create a df with connectance, extinctions and variance of extinctions
data.ext_con <- con     #set up df, first col will be connectance
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, mean))) #add mean no. of extinctions of n=rep iterations
data.ext_con <- rbind(data.ext_con, t(apply(extinctions_mat, MARGIN = 2, var)))  #add variance of extinctions of n=rep iterations
data.ext_con <- data.ext_con %>% #transpose and save as dataframe
  t() %>%
  as.data.frame()
  
#rename rows and columns:
rownames(data.ext_con) <- NULL
colnames(data.ext_con) <- c("con", "ext.m", "ext.v")

#plot the mean extinction numbers
ggplot(data.ext_con, aes(x = con, y=ext.m))+
  geom_point()+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.m)+1)) #becomes asymptotic to x-axis parallel with connectance -> 0.5

#plot the variance of the extinctions
ggplot(data.ext_con, aes(x = con, y=ext.v))+
  geom_point()+
  scale_y_continuous(limits = c(0, max(data.ext_con$ext.v)+1)) #interesting pattern, lets have a look at it with 10000 rep's


####4.2 plotting the same stuff for a subset food web####
NULL






