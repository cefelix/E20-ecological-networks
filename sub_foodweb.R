####what do i want to do?####

#explanatory1: temp
#explanatory2: connectance

#response1: individual abundances /biomass

#select subset species based on:
  #different number:
  #n_bas = c(1,2,4,8,16) 


#Q default is 1.2, set it to 1
#model$C <- 0
#don't go for significance
#connectance  from 0.5 to 0.4 in 0.05 steps



#libraries:
library(ATNr)
library(dplyr)

####1.1 exp = temp; res = n_extinct####

n_bas = 8
n_tot = 32 #use a multiple of 8
n_sub = seq(0,1, by=0.125)*(n_tot-n_bas) 
  #set subset_foodweb to 25% (respectively: 50,75,100 %) of total species, select randomly


#create an ordered body mass vector:
set.seed(666)
BM <- runif(n_tot, 2, 3) %>%
  sort()
BM <- 10^BM
plot(BM)

#L matrix
L <- create_Lmatrix(BM, n_bas, Ropt = 50) # for


#create a food web matrix that contains 0/1
fw.1 <- L
fw.1[fw.1>0.0] <- 1

#fw.1 <- create_niche_model(n_tot, 0.4)

####1.2 model creation#####
model_t.ex <- create_model_Unscaled(n_tot, n_bas, BM, fw.1) %>%
  initialise_default_Unscaled(temperature = 40)
#model_t.ex$c <- rep(0, (n_tot-n_bas))
#model_t.ex$q <- rep(1, (n_tot-n_bas))
model_t.ex$q %>% summary
model_t.ex$c %>% summary



#biomass dynamics:
times <- seq(1, 1e8, by = 1e6) #time for integration
biomass <- runif(n_tot, 1, 4) #initial biomasses for each species
sol = lsoda_wrapper(times, biomass, model = model_t.ex)
sol %>%
  str()

plot_odeweb(sol, n_tot) #seems as there are almost no trophic interactions



#looking at extinctions
finals = sol[nrow(sol), -1]
sum(finals <= model_t.ex$ext)


####1.3 iterating multiple times over a temperature gradient####
 
#times <- seq(1, 1e4, by = 1e2) #this timeseries runs too short! try 1e8!
temperatures <- seq(10, 30, by= 2)
i = 0
extinctions = rep(NA, length(temperatures))
iterations = 20

extinctions_out <- c()
temperatures_out <- c()

biomasses_out_0.5 <- c() #lets check the ratios of biomasses at the middle/end of time
                         #(this should become a list where i can append the vectors to)
biomasses_out_1 <- c()   #as above


for (i in 1:iterations) {
  for (temp in temperatures) {
    print(temp)
    model = initialise_default_Unscaled(model_t.ex, temperature = temp)
    sol=lsoda_wrapper(times, biomass, model)
    finals = sol[nrow(sol), -1]
    exts = sum(finals < model$ext)
    
    biom_0.5 = sol[nrow(sol)/2, -1]
    biomasses_out_0.5 = c(biomasses_out_0.5, biom_0.5)
    biom_1 = sol[nrow(sol), -1]
    biomasses_out_1 = c(biomasses_out_1, biom_1)
    
    
    extinctions_out = c(extinctions_out, exts)
    temperatures_out= c(temperatures_out, temp)
    
  }
  i = i+1
}

plot(extinctions_out ~ temperatures_out)  #no extinctions o.O -> model doesn't work?

sum((biomasses_out_0.5 == biomasses_out_1)==TRUE) # so there is some dynamics
sum((abs(biomasses_out_0.5 - biomasses_out_1) > 0.1) == TRUE)


####2.1 explanatory: Ropts using an allometric model####

#explanatory:
  Ropt <- seq(10, 70, by=15)

#responses (each to be compared, between a food web of n species and a sub-foodweb of k<n species)
  No.ext #number of extinctions
  abundances #total biomass / Body mass
  connectance #in the final food web
  
  
  
#iterate through 20 replicates per food web  

####3.1 explanatory: connectance using a niche-model####
  
####3.2 defining explanatory/response variables####  
  #explanatory:
  connectance <- seq(0.2, 0.45, by=0.05)
  
  #responses (each to be compared, between a food web of n species and a sub-foodweb of k<n species)
  No.ext #number of extinctions
  abundances #total biomass / Body mass
  connectance #in the final food web

####3.3 set up parameters####
  
  