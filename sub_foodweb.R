####what do i want to do?####

#explanatory1: temp
#explanatory2: connectance

#response1: n.extinctions, individual abundances /biomass

#select subset species based on:
  #different number:
  #n_bas = c(1,2,4,8,16) 


#libraries:
library(ATNr)


####1.1 exp = temp; res = n_extinct####

n_bas = 8
n_tot = 32 #use a multiple of 8
n_sub = seq(0,1, by=0.125)*(n_tot-n_bas) 
  #set subset_foodweb to 25% (respectively: 50,75,100 %) of total species, select randomly


#create an ordered body mass vector:
set.seed(123)
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
  initialise_default_Unscaled(mod)


#biomass dynamics:
times <- seq(1, 1e3, 1e2) #time for integration
biomass <- runif(n_tot, 2, 6) #initial biomasses for each species
sol = lsoda_wrapper(times, biomass, model = temp_ext)
sol %>%
  str()

plot_odeweb(sol, n_tot) #seems as there are almost no trophic interactions



