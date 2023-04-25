#statistical analysis
library(ATNr)

data <- read.csv("./data/128spec.csv")[,-1] #drops first column which is rownumber
d <- data

d_real <- subset(d, d$con <= 0.30)
d <- d_real
####scatterplots extinctions####


plot(d$con, d$ext_all)
  plot(d$ext_big, d$ext_all)
  plot(d$ext_small, d$ext_all)
plot(d$con, d$ext_big)
plot(d$con, d$ext_small)

####scatterplots shannon####
plot(d$con, d$shan_all)
  plot(d$shan_big, d$shan_all)
  plot(d$shan_small, d$shan_all)
plot(d$con, d$shan_big)

####scatterplots trophic level####
plot(d$con, d$minTL_B) 
summary(d$minTL_s)
plot(d$con, d$shan_small)

log(32)#
log(128) #maximum shannon index 

