####what do i want to do?####

#explanatory1: temp
#explanatory2: connectance

#response1: n.extinctions, individual abundances /biomass

#select subset species based on:
  #different number:
  #n_bas = c(1,2,4,8,16) 


#libraries:
library(ATNr)
library(dplyr)

####1.1 exp = temp; res = n_extinct####

n_bas = 8
n_tot = 100 #use a multiple of 8
n_sub = seq(0,1, by=0.125)*(n_tot-n_bas) 
  #set subset_foodweb to 25% (respectively: 50,75,100 %) of total species, select randomly


#create an ordered body mass vector:
set.seed(187023)
BM <- runif(n_tot, 2, 9) %>%
  sort()
BM <- 10^BM
plot(BM)

#L matrix
L <- create_Lmatrix(BM, n_bas)


#create a food web matrix that contains 0/1
fw.1 <- L
fw.1[fw.1>0.0] <- 1

####1.2 model creation#####
temp_ext <- create_model_Unscaled(n_tot, n_bas, BM, fw.1) %>%
  initialise_default_Unscaled(temperature = 40)


#biomass dynamics:
times <- seq(1, 1e4, by = 1e2) #time for integration
biomass <- runif(n_tot, 2, 6) #initial biomasses for each species
sol = lsoda_wrapper(times, biomass, model = temp_ext)
sol %>%
  str()

plot_odeweb(sol, n_tot) #seems as there are almost no trophic interactions

finals = sol[nrow(sol), -1]
sum(finals <= temp_ext$ext)


####1.3 iterating multiple times over a temperature gradient####

times <- seq(1, 1e4, by = 1e2)
temperatures <- seq(10, 30, by= 2)
i = 0
extinctions = rep(NA, length(temperatures))
iterations = 20

extinctions_out <- c()
temperatures_out <- c()

for (i in 1:iterations) {
  for (temp in temperatures) {
    print(temp)
    model = initialise_default_Unscaled(temp_ext, temperature = temp)
    sol=lsoda_wrapper(times, biomass, model)
    finals = sol[nrow(sol), -1]
    exts = sum(finals < model$ext)
    extinctions_out = c(extinctions_out, exts)
    temperatures_out= c(temperatures_out, temp)
    
  }
  i = i+1
}

plot(extinctions_out ~ temperatures_out)  #no extinctions o.O -> model doesn't work
