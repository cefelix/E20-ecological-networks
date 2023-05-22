#statistical analysis
library(ATNr)
library(dplyr)
library(vegan)
library(ggplot2)


d <- read.csv("./data/20220506_96spec_2000c_04to12bas_v02.csv")[,-1] #drops first column which is rownumber

#### 0 - general influence of connectance on extinctions and diversity####
plot(d$con, d$ext_all)

model.all <- lm(ext_all ~ con, data = d)
  summary(model)
exts_pred <- data.frame(pred.exts = predict(model, d)) 
  
#fig 1a  
ggplot(data=d, aes(x=con, y=ext_all))+
  geom_point()+
  #geom_line(color="red", data = exts_pred, aes(x=con, y=pred.exts))+
  theme_classic()+
  scale_y_continuous(name= "Extinctions")+
  scale_x_continuous(name= "Connectance", breaks = seq(0.05, 0.35, by=0.05), limits=c(0.05, 0.35))


#fig 1b
ggplot(data=d, aes(x=con, y=shan_all))+
  geom_point()+
  #geom_line(color="red", data = exts_pred, aes(x=con, y=pred.exts))+
  theme_classic()+
  scale_y_continuous(name= "Shannon Index")+
  scale_x_continuous(name= "Connectance", breaks = seq(0.05, 0.35, by=0.05), limits=c(0.05, 0.35))

#fig 1c
#using shannon equitability index instead
d$shan_equ <- d$shan_all / log(96-d$ext_all)
ggplot(data = d, aes(x=con, y=shan_equ))+
  geom_point()+
  theme_classic()+
  scale_y_continuous(name= "Shannon Equitability Index")+
  scale_x_continuous(name= "Connectance", breaks = seq(0.05, 0.35, by=0.05), limits=c(0.05, 0.35))



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
  
  sum(extinctions[rowSums(extinctions) > 96-24]) # 1589 -> only ~400 fw's with less th
  
#producing nice plots now
ggplot(d, aes(ext_all, ext_big))+
  geom_point()
  
    
ggplot(d, aes(shan_all, shan_big))+
  geom_point()
  
  
#### 2 - comparing change of correlation in the randomly sampled species along the gradient in sampling size####

####2.1 data preparation####  
#extinctions  
d_cor <- c(cor(d$ext_all, d$rand_ext_4),
           cor(d$ext_all, d$rand_ext_8),
           cor(d$ext_all, d$rand_ext_12),
           cor(d$ext_all, d$rand_ext_16),
           cor(d$ext_all, d$rand_ext_20),
           cor(d$ext_all, d$rand_ext_24),
           cor(d$ext_all, d$rand_ext_28),
           cor(d$ext_all, d$rand_ext_32),
           cor(d$ext_all, d$rand_ext_36) 
          )

d_cor <- cbind(d_cor,
               c(cor(d$ext_all, d$rand_extTL_4),
                 cor(d$ext_all, d$rand_extTL_8),
                 cor(d$ext_all, d$rand_extTL_12),
                 cor(d$ext_all, d$rand_extTL_16),
                 cor(d$ext_all, d$rand_extTL_20),
                 cor(d$ext_all, d$rand_extTL_24),
                 cor(d$ext_all, d$rand_extTL_28),
                 cor(d$ext_all, d$rand_extTL_32),
                 cor(d$ext_all, d$rand_extTL_36) 
                 )
                )   


  
#shannon indices  
d_cor <- cbind(d_cor,
               c(cor(d$shan_all, d$rand_shan_4),
                 cor(d$shan_all, d$rand_shan_8),
                 cor(d$shan_all, d$rand_shan_12),
                 cor(d$shan_all, d$rand_shan_16),
                 cor(d$shan_all, d$rand_shan_20),
                 cor(d$shan_all, d$rand_shan_24),
                 cor(d$shan_all, d$rand_shan_28),
                 cor(d$shan_all, d$rand_shan_32),
                 cor(d$shan_all, d$rand_shan_36) 
                 )
                ) 

d_cor <- cbind(d_cor,
               c(cor(d$shan_all, d$rand_shanTL_4),
                 cor(d$shan_all, d$rand_shanTL_8),
                 cor(d$shan_all, d$rand_shanTL_12),
                 cor(d$shan_all, d$rand_shanTL_16),
                 cor(d$shan_all, d$rand_shanTL_20),
                 cor(d$shan_all, d$rand_shanTL_24),
                 cor(d$shan_all, d$rand_shanTL_28),
                 cor(d$shan_all, d$rand_shanTL_32),
                 cor(d$shan_all, d$rand_shanTL_36)
                 )
                )

#colnames, adding sample size

d_cor <- cbind(seq(4, 36, by=4), d_cor)
  d_cor <- as.data.frame(d_cor)
  colnames(d_cor) <- c("sample_size" ,"cor_ext_rand", "cor_ext_randTL", 
                     "cor_shan_rand", "cor_shan_randTL")
  
write.csv(d_cor, "./data/20220517_correlations_v01.csv")  
  
d_cor <- read.csv("./data/20220517_correlations_v01.csv")  


#plotting correlation coefficients  
ggplot(data=d_cor, aes())+
  geom_point(aes(x = sample_size, y=cor_shan_rand, col= "red"))+
  geom_point(aes(x = sample_size, y=cor_shan_randTL, col= "pink"))+
  geom_point(aes(x = sample_size, y=cor_ext_rand, col= "blue"))+
  geom_point(aes(x = sample_size, y=cor_ext_randTL, col= "black"))+
  scale_y_continuous(name = "pearson correlation coefficient",limits = c(0,1))+
  scale_x_continuous(limits = c(4,36), breaks= seq(4,36, by=8))+
  theme_classic()
  




####2.2 analysis####
  plot(d_cor$sample_size, d_cor$cor_ext_rand)
  plot(d_cor$sample_size, d_cor$cor_ext_randTL)

  plot(d_cor$sample_size, d_cor$cor_shan_rand) #fluctuations at lower sample sizes due to calculation of shannon index including extinct species
  plot(d_cor$sample_size, d_cor$cor_shan_randTL)
  
ggplot(data=d_cor, )
  
  
####2.3 plotting####
  knitr::kable(d_cor)


####3.1 comparing big, small, random, random with TL variation samples of 24 species####
d <- read.csv("./data/20230522_SampleSizes_v01.csv")[,-1]
d$subsize <- factor(d$subsize, levels=c("4", "8", "12", "16", "20", 
                                        "24", "28", "32", "36")  )
d <- subset(d, subsize == "24")

n_cons <- rep(NA, length(con))
for (i in 1:length(con)) {
  n_bas <- sum(troph.lvl[i]== 1)
  n_cons[i] <- 96-n_bas
  
}
d$n_cons <- n_cons
#adding extinction rate ratios in big /small consumer samples
d$ext_ratio_sml <- (extinctions.small_vec / 24 ) / (extinctions_vec / d$n_cons)
d$ext_ratio_big <- (extinctions.big_vec / 24)    / (extinctions_vec / d$n_cons)
                    #extinction rate in the sample / #extinction rate in whole food web

n_cons <- rep(NA, length(con))
for (i in 1:length(con)) {
  n_bas <- sum(troph.lvl[i]== 1)
  n_cons[i] <- 96-n_bas
  
}
d$n_cons <- n_cons

d.long <- pivot_longer(d, cols = c("ext_ratio", "ext_ratioTL", "ext_ratio_sml", "ext_ratio_big"))
str(d.long)
colnames(d.long) <- c("con", "subsize", "n_consumers", "sampling", "ratio")
d.long$sampling <- factor(d.long$sampling, 
                          levels = c("ext_ratio_sml", "ext_ratio_big",
                                     "ext_ratio", "ext_ratioTL"))
levels(d.long$sampling) <- c("lowest BM", "highest BM", "random", "TroLev")
str(d.long)
  #done with data arranging, now plotting:

write.csv(d.long, "./data/20230523_fig2a_24BiSmRaTL.csv")
d.long <- read.csv("./data/20230523_fig2a_24BiSmRaTL.csv")

ggplot(d.long)+
  geom_jitter(aes(y=ratio, x=sampling), size=0.5, alpha=0.3)+
  stat_summary(aes(y=ratio, x=sampling), fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", color="red", width=0.3, linewidth=0.8)+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  ggtitle("Comparison of sampling strategies for 24 species")+
  scale_y_continuous(name = "Extinction-rate ratio", limits = c(-0.2, 8))+
  scale_x_discrete("Sampling strategy")+
  theme_classic()

subset(d.long, sampling== "random")$ratio %>%
  sd() #mean=1.083, sd=0.336
subset(d.long, sampling== "TroLev")$ratio %>%
  sd() #mean=1.093, sd=0.351 

subset(d.long, sampling== "lowest BM")$ratio %>%
  mean() #mean=0.621, sd=0.446
subset(d.long, sampling== "highest BM")$ratio %>%
  sd()   #mean=1.708, sd=0.538



mean(d.long$ratio)


####3.2 jitter plots of variation in extinction rate between sample sizes####
d <- read.csv("./data/20230522_SampleSizes_v01.csv")[,-1]
d$subsize <- factor(d$subsize, levels=c("4", "8", "12", "16", "20", 
                                        "24", "28", "32", "36")  )
str(d)

summary(is.na(d$ext_ratioTL)) #no NAs
summary(is.na(d$ext_ratio))   #no NAs



#random sampling
ggplot(d) +
  geom_jitter(aes(y=ext_ratio, x=subsize), size=0.5)+
  stat_summary(aes(y=ext_ratio, x=subsize), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.3, linewidth=0.8)+
  stat_summary(aes(y=ext_ratio, x=subsize), fun = mean, geom = "point",
               color="red")+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  scale_x_discrete(name = "No. of consumer species sampled")+
  scale_y_continuous(name="Extinction-rate ratio", limits=c(-0.2,8))+
  ggtitle("Sampling randomly (fig. 2a)")+
  theme_classic()


#random sampling spanning several trophic levels
ggplot(d) +
  geom_jitter(aes(y = ext_ratioTL,   x=subsize), size=0.5, alpha=0.5)+
  stat_summary(aes(y=ext_ratioTL, x=subsize), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width= 0.3, linewidth=0.8)+
  stat_summary(aes(y=ext_ratioTL, x=subsize), fun = mean, geom = "point",
               color="red")+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  scale_x_discrete(name = "No. of consumer species sampled")+
  scale_y_continuous(name="Extinction-rate ratio", limits=c(-0.2,8))+
  ggtitle("Sampling from different trophic levels (fig. 2b)")+
  theme_classic()


#compare random sampling and trophic level dependent sampling at 4, 20 and 36 sampled species
d.long <- pivot_longer(d, cols = c("ext_ratio", "ext_ratioTL"))
colnames(d.long) <- c("con", "subsize", "sampling", "ratio")
d.long$sampling <- factor(d.long$sampling)
levels(d.long$sampling) <- c("random", "TroLev")
str(d.long)

d.long <- subset(d.long, subsize %in% c("4", "20", "36"))

write.csv(d.long, "./data/20230522_dat_fig3_v01.csv")
d.long <- read.csv("./data/20230522_dat_fig3_v01.csv")

ggplot(d.long)+
  geom_jitter(aes(y=ratio, x=sampling), size=0.5, alpha=0.3)+
  stat_summary(aes(y=ratio, x=sampling), fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", color="red", width=0.3, linewidth=0.8)+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  facet_wrap(vars(subsize), labeller = labeller( subsize=
                                                   c("4" = "4 species", "20" = "20 species", 
                                                     "36" = "36 species")
                                                 )   )+
  ggtitle("Comparison of sampling strategies")+
  scale_y_continuous(name = "Extinction-rate ratio", limits = c(-0.2, 8))+
  scale_x_discrete("Sampling strategy at different sample sizes")+
  theme_classic()


####3.3 jitter plots if variation in shannon index at different sample sizes####
d <- read.csv("./data/20230523_SampleSizes_v01.csv")[,-1]
d$subsize <- factor(d$subsize, levels=c("4", "8", "12", "16", "20", 
                                        "24", "28", "32", "36")  )
str(d)

#random sampling
ggplot(d) +
  geom_jitter(aes(y=shan_TL, x=subsize), size=0.5)+
  stat_summary(aes(y=shan_TL, x=subsize), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.3, linewidth=0.8)+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  scale_x_discrete(name = "No. of consumer species sampled")+
  scale_y_continuous(name="Shannon equitability index ratio")+
  ggtitle("Sampling randomly (fig. 2a)")+
  theme_classic()

#TL sampling
ggplot(d) +
  geom_jitter(aes(y=shan_rand, x=subsize), size=0.5)+
  stat_summary(aes(y=shan_rand, x=subsize), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.3, linewidth=0.8)+
  geom_hline(yintercept = 1, linetype="dashed", color="cadetblue", linewidth=1.0)+
  scale_x_discrete(name = "No. of consumer species sampled")+
  scale_y_continuous(name="Shannon equitability index ratio")+
  ggtitle("Sampling randomly (fig. 2a)")+
  theme_classic()
