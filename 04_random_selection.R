#select randomly but based on conditions
library(ATNr)


#read data (output of code #3)
d <- read.csv("./data/64spec_2000randoms.csv")[,-1] #drops first column which is rownumber

#read RDS files (output of code #2): (DELETE ?! - irrelevant)
abundance_array <- readRDS(file = "./raw/abundance.rds")
biomass_array <- readRDS(file = "./raw/biomass.rds")
extinction_array <-  readRDS(file = "./raw/extinction.rds")
troph.lvl_array <- readRDS(file = "./raw/trophic_lvl.rds")


####
str(d)

####extinctions####
#condition = 

while (condition) {
  
}

