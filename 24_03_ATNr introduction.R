library(devtools)
install_github('gauzens/ATNr')

library(ATNr)
library(dplyr)
set.seed(123)


####define community of species####
nb_s = 30 #number of species
nb_b = 8 #basal species
nb_n = 2 #number of nutrient pools
  
BM = runif(30, 2, 9)
BM = sort(BM)
BM = 10^BM
BM

####define trophic interactions####

L = create_Lmatrix(BM, nb_b)
  #Ropt = peak of distribution in nature fig 2
  #gamma defines width of this distribution
L %>% 
  str()


#create a food web matrix that contains 0/1
fw = L
fw[fw>0.0] = 1
fw

model <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, BM, fw)
model %>% 
  str()
model$BM #the same as BM

model = initialise_default_Unscaled_nuts(model, L)
model %>%
  str()


times = seq(1, 1e4, 1e2) #time for integration
bioms= runif(nb_s + nb_n, 2, 3) #intial biomasses and nutrient concentrations
sol = lsoda_wrapper(times, bioms, model)
sol %>%
  str()

plot_odeweb(sol, nb_s) 

#look at possible extinctions --> look at final densities, last row of sol df 
#remove the collumns with time and nutrients
sol[nrow(sol),-c(1,2,3)]
#or, more flexible:

finals = sol[nrow(sol), -(1:(nb_n+1))]

#count number of non-extinct species:
sum(finals > 0)
#or better use model treshold below which species would be considered extinct
sum(finals < model$ext)

#number of ext basal species
sum(finals[1:nb_b] < model$ext)

#number of ext non-basal species
sum(finals[(nb_b+1):nb_s] < model$ext)

#https://cran.r-project.org/web/packages/ATNr/index.html

####redo it for n models####

nb_s = 30 #number of species
nb_b = 8 #basal species
nb_n = 2 #number of nutrient pools

BM = runif(30, 2, 9)
BM = sort(BM)
BM = 10^BM
BM

####define trophic interactions####

L = create_Lmatrix(BM, nb_b)
#Ropt = peak of distribution in nature fig 2
#gamma defines width of this distribution
L %>% 
  str()


#create a food web matrix that contains 0/1
fw = L
fw[fw>0.0] = 1
fw

model <- create_model_Unscaled_nuts(nb_s, nb_b, nb_n, BM, fw)
model %>% 
  str()
model$BM #the same as BM

#look over a temperature gradient
temperatures <- seq(10,30, by = 2)
times = seq(1, 1e4, 1e2) #time for integration
bioms= runif(nb_s + nb_n, 2, 3) #intial biomasses and nutrient concentrations

extinctions = rep(NA, length(temperatures))
i=0

replicates = 10

extinctions = c()
temps = c()

for(rep in 1:replicates){
  for(temperature in temperatures){
    i = i+1
    print(temperature)
    model = initialise_default_Unscaled_nuts(model, L, temperature = temperature)
    sol=lsoda_wrapper(times, bioms, model)
    finals = sol[nrow(sol), -(1:(nb_n+1))]
    exts = sum(finals < model$ext)
    extinctions = c(extinctions, exts)
    temps= c(temps, temperature)
  }
} #using a distrubution of extinctions instead of a point aestimate
#if change parameters inside model, recreate model inside loop!!!



#for(temperature in temperatures){
  i = i+1
  print(temperature)
  model = initialise_default_Unscaled_nuts(model, L, temperature = temperature)
  sol=lsoda_wrapper(times, bioms, model)
  finals = sol[nrow(sol), -(1:(nb_n+1))]
  exts = sum(finals < model$ext)
  extinctions[i] = exts
} #extinctions are point aestimates


plot(extinctions ~ temps)

model = initialise_default_Unscaled_nuts(model, L, temperature = temperatures)


times = seq(1, 1e4, 1e2) #time for integration
bioms= runif(nb_s + nb_n, 2, 3) #intial biomasses and nutrient concentrations
sol = lsoda_wrapper(times, bioms, model)
sol %>%
  str()

#####alternative
?create_niche_model


c <- seq(1,4)
a <- matrix(c,2,2)

a[1,] <- c(8,6)
a[,1] <- c(8,6)
