#statistical analysis
library(ATNr)

d <- read.csv("./data/20220506_96spec_2000c_04to12bas_v02.csv")[,-1] #drops first column which is rownumber

#### 1 - comparing pearson correlations in 24 species samples####


#extinctions
cor(d$ext_all, d$ext_big)           #0.868
cor(d$ext_all, d$ext_small)         #0.454
cor(d$ext_all, d$rand_ext_24)       #0.912
cor(d$ext_all, d$rand_extTL_24)     #0.923

cor(d$ext_all ,d$ext_BAS) #sd = zero
  sum(d$ext_BAS)          #0

  
#diversity
cor(d$shan_all, d$shan_big)         #0.427
  plot(d$shan_all, d$shan_big)
  sum(d$shan_all < 0.25)
  sum(d$shan_big < 0.25)
  
cor(d$shan_all, d$shan_small)       #0.390
  plot(d$shan_all, d$shan_small)
  sum(d$shan_small < 0.25)

cor(d$shan_all, d$rand_shan_24)     #0.249
  plot(d$shan_all, d$rand_shan_24)
  sum(d$rand_shan_24 < 0.25)

cor(d$shan_all, d$rand_shanTL_24)   #0.299
  plot(d$shan_all, d$rand_shanTL_24)
  sum(d$rand_shanTL_24 < 0.25)
  
  
#### 2 - comparing change of correlation in the randomly sampled species along the gradient in sampling size####
  
  
  