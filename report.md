# Title:How does a selection of species represent parameters of the entire food web based on selection criteria?

## Abstract

## Introduction

Many ecological parameters are not measured directly, but assessed through indicator species.

Food webs blabla ...

## 2 Methods

I did all data generation and analysis using `R 4.2.2` (R Core Team, 2022). Food web dynamics were simulated using the `ATNr` package by Gauzens et al. (2022). Other packages used are `vegan` (Oksanen et al., 2022), `ggplot2` (Wickham, 2016), and `tidyr` (Wickham et al., 2023). The complete code for this project can be found at [github.com/cefelix/e20-ecological-networks](https://github.com/cefelix/E20-ecological-networks).

### 2.1 Simulation of food web dynamics

I simulated food web dynamics for 2000 food webs, each food web consisting of `n_tot = 96` trophic species. Connectance `con` varied from 0.05 to 0.35, which corresponds roughly to the span of connectance in empirical food webs as reported by Dunne et al. (2002). Also, I set up the variables `times`, `biom`, and `ext_tresh` which are explained later in this chapter. Finally, one array per output variable was set up using `array(dim = c(length(con), n_tot))`. Exemplary I show the array to store final biomasses below.

``` r
n_tot = 96                                          #total number of species
con <- runif(2000, min=0.05, max=0.35)  %>% sort()  #2000 connectance values

times <- seq(1, 1e12, by = 1e9)   #time for integration of biomass dynamics
biom <- runif(n_tot, 1, 4)        #initial biomasses
ext_thresh <- 0.1**6              #threshold below which considered extinct

biomasses <- array(dim = c(length(con), n_tot)) #empty array
```

Next, I generated an adjacency matrix `fw` generated using `ATNr::create_niche_model()`. I assessed the number of basal species `n_bas` by calculating `colSums(fw)`. To minimize effects of variation in basal species number and thus disentangle effects of variation in connectance from effects of variation in basal species number, adjacency matrices were restricted to have between 4 and 12 basal species. This arbitrary range helped avoiding extreme numbers of basal species (the range without restriction lay between 1 and 20+ basal species, with food webs of low connectance typically having very few basal species, while highly connected ones having plenty). Then I set up a body mass vector `BM` with each element representing the body mass of one species. `BM` spanned 12 orders of magnitude, which is considered realistic in soil food webs (Potapov, 2021). All following steps in chapters 2.1 and 2.2 are to be carried out within the loop through `con` (shown in the code section below).

``` r
for(i in 1:length(con)) {
  
  fw <- create_niche_model(S = n_tot, C = con[i])
  n_bas <- sum(colSums(fw) == 0)
    
  while (n_bas > 12 | n_bas < 4) { #otherwise create new adjacency matrix
    fw <- create_niche_model(S = n_tot, C = con[i])
    n_bas <- sum(colSums(fw) == 0)
  }
  
  BM <- runif(n_tot, 1, 12) %>% 
      sort()
  BM <- (10^BM) 
  #[...]
}
```

Next, an allometric trophic model `model` was created using `create_model_Unscaled()`. Parameters for the ordinary differential equation (ODE) system describing biomass dynamics for each species were initialised by `initialise_model_Unscaled()`. The extinction threshold `model$ext` (i.e. the biomass below which a species is considered extinct) was set to 1e-6 manually. ODE's were then solved numerically over a pre-defined time course `times` using `lsoda_wrapper()`, and the final biomasses were stored in `sol`.

``` r
model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
  initialise_default_Unscaled()
model$ext <- ext_thresh

sol <- lsoda_wrapper(times, biom, model) 
```

### 2.2 Extraction of output variables

I extracted the final biomasses `bioms` for each species using `sol[nrow(sol), -1]` and stored them in `biomasses[i, ]` (the first dimension of `biomasses` corresponds to connectance, while the second corresponds to the species). I applied the same principle to store all the other output variables.

``` r
bioms = sol[nrow(sol),-1]       #final biomasses for each species
biomasses[i,] = bioms           #storing in array at i-th element of connectance
```

Abundances were calculated dividing `bioms` by individual body masses: `abuns = bioms / BM`. Extinction of each single species was assessed by comparing `bioms` to the previously set extinction threshold 1e-6 using `exts = ifelse(bioms <= model$ext, yes = 1, no = 0)`. Trophic level for each species was calculated with `TroLev(fw)` according to the general formula *INSERT*. The number of initial trophic links (number of prey species `prey.on_start` and consumer species `consumed.by_start`) for each species was calculated by `colSums(fw)` and `rowsums(fw)`, respectively.

``` r
abuns = bioms / BM                                  #final abundances
exts = ifelse(bioms <= model$ext, yes = 1, no = 0)  #extinctions
trolev =  TroLev(fw)                                #trophic levels
prey.on_start = colSums(fw)                         #trophic links: predator
consumed.by_start = rowSums(fw)                     #trophic links: consumed 
```

I then cleared the initial adjacency matrix `fw` from the interactions that involved at least one extinct species. Therefore, I set rows and columns of `fw` to zero in case the corresponding species had become extinct, and stored the result in `fw_end`. The number of final trophic links were calculated as `prey.on_end = colSums(fw_end)` and `consumed.by_end = rowSums(fw_end)`. Afterwards, I stored all output variables in the corresponding arrays, resulting each array to store parameters for 96 species at 2000 unique values of connectance.

``` r
fw_end <- fw
fw_end[exts == 1,] <- 0       #if extinct species: row is set to zero
fw_end[,exts == 1] <- 0       #if extinct species: column is set to zero

prey.on_end = colSums(fw_end)               #trophic links: predator          
consumed.by_end = rowSums(fw_end)           #trophic links: consumed

#storing output in arrays:

biomasses[i,]         = bioms 
abundances[i,]        = abuns 
extinctions[i,]       = exts
troph.lvl[i,]         = trolev

prey_on_START[i,]     = prey_links 
consumed_by_START[i,] = pred_links 
prey_on_END[i,]       = prey.on_start
consumed_by_END[i,]   = consumed.by_start
```

### 2.3 Sampling from sub-food webs: 24 smallest / 24 biggest consumers

I calculated extinctions in each single food web by applying `sum()` over each row of the `extinctions` array and saving the result in `extinctions_vec`. Analogously, I calculated Shannon indices for each food web by applying `vegan::diversity()` over each row the `abundances` array [^report-1] and saved the result in `shannon_vec`. Furthermore, I set up empty vectors to store extinctions in the basal species `extinctions.BAS_vec`, the smallest consumer species `extinctions.small_vec`, and the biggest consumer species `extinctions.big_vec`. Analogously, I set up empty vectors to store Shannon indices. Lastly, I set up another empty vector `n.bas_vec` to store the number of basal species in each food web.

[^report-1]: Technically, the here calculated abundances are densities. It is debated whether it is appropriate to calculate Shannon-indices from densities, but as this is not the focus of this report, I just went ahead with it.

``` r
extinctions_vec <- apply(extinctions, MARGIN = c(1), FUN = sum) %>% 
  as.vector()       #extinctions in each entire food web
  extinctions.BAS_vec     <- rep(NA, length(con))
  extinctions.small_vec   <- rep(NA, length(con))
  extinctions.big_vec     <- rep(NA, length(con))

shannon_vec <-  apply(abundances, MARGIN = c(1), FUN = diversity) %>%
  as.vector()       #shannon index in each entire food web
  shannon.BAS_vec     <- rep(NA, length(con))  
  shannon.small_vec   <- rep(NA, length(con))
  shannon.big_vec     <- rep(NA, length(con))
  
n.bas_vec <- rep(NA, length(con))  
```

In order to select the 24 biggest and 24 smallest consumer species, I initialized two selection vectors `slct_bi` and `slct_sm`. This is straightforward, as species are already sorted by bodymass when feeding the sorted vector `BM` into `create_model_Unscaled(n_tot, n_bas, BM, fw)`. As the number of basal species varied across food webs, this had to be done for each food web separately using another for-loop. As the first dimension of `extinctions` indicates the value of connectance `con` (which corresponds to the first simulated food web) and the second dimension indicates the species, calling `extinctions[j, slct_sm]` prints a binary vector, showing whether each of the smallest consumer species went extinct (1) or not (0). Applying `sum()` onto this vector, gives the total number of extinctions within the selection of the smallest consumer species `slct_sm`. This is stored at the j-th position in the previously set up `extinctions.small_vec`. Similarily, `diversity()` is applied on the subset of the `abundances` array.

``` r
n_tot = 96                                    #total number of species per food web
n_sub = 24                                    #no. of species in slct_sm / slct_bi

for(j in 1:length(con)) {
    
    bas <- n_tot - sum(troph.lvl[j,] > 1)     # number of basal species
    consumers <- c((bas+1):n_tot)             # a vector of all consumer species
    n.bas_vec[j] <- bas                       # number of basal species 
    
    slct_sm <- consumers[1:n_sub]             # the smallest consumer species
    slct_bi <- c((n_tot-n_sub+1):n_tot)       # the biggest consumer species

    #extinctions in 24 species
    extinctions.small_vec[j] <- sum(extinctions[j, slct_sm]) #smallest consumers
    extinctions.BAS_vec[j] <- sum(extinctions[j, c(1:bas)])  #basal species
    extinctions.big_vec[j] <- sum(extinctions[j, slct_bi])   #biggest consumers

    #diversity in 24 species
    shannon.small_vec[j] <- diversity(abundances[j, slct_sm])   #smallest consumers
    shannon.BAS_vec[j]   <- diversity(abundances[j, c(1:bas)])  #basal species
    shannon.big_vec[j] <- diversity(abundances[j, slct_bi])     #biggest consumers
}
```

### 2.5 Sampling from sub-food webs: random species selection

I also selected from 4 up to 36 random species from each food web. For this, I initialized a vector `selection_size <- seq(from = 36, to = 4, by=-4)` ranging from 4 to 36 by steps of 4, indicating the number of random species I would select. I then set up empty arrays to store the number of extinctions and shannon indices within the random selection of species for each sample size.

``` r
selection_size          <- seq(from = 36, to = 4, by=-4) 
extinctions.rand_mat    <- array(NA, dim = c( length(con), length(selection_size) ))
shannon.rand_mat        <- array(NA, dim = c( length(con), length(selection_size) )) 
```

I then looped through all food webs using the outer for-loop `for (j  in 1:length(con))`. For each food web, I calculated a vector containing all consumer species `consumers` the same way as in chapter 2.4. I then looped through the `selection_size` vector and sampled a random selection of consumer species using `sample()` with `size = selection size`. The number of extinctions and Shannon indices within this selection were calculated analogously to chapter 2.4 and stored in extinctions.rand_mat and shannon.rand_mat, respectively.

``` r
for (j  in 1:length(con)) {
  bas <- n_tot - sum(troph.lvl[j,] > 1)   # number of basal species
  consumers <- c((bas+1):n_tot)           # a vector of all consumer species
  
  for (k in 1:length(selection_size)) {
    slct_rand <- sample(consumers, size = selection_size[k]) %>% 
      sort()                              #a random selection of 4 to 36 species
    
    extinctions.rand_mat[j,k] <- sum(extinctions[j, slct_rand])     
    shannon.rand_mat[j,k] <- diversity((abundances[j, slct_rand]))  
  }
}
```

### 2.6 Sampling from sub-food webs: random species selection, but forcing variation in trophic levels

Lastly, I sampled 4 to 36 random consumer species from each food web, but this time restricting the random selection to span several trophic levels. To achieve this, `slct_rand` was set up differently: For each food web, I extracted the trophic levels of the consumer species by first accessing `troph.lvl[j,consumers]` and then storing the output in `TL_c`. Afterwards I sampled randomly from each quarter of consumer species (sorted by their trophic levels): I compared `TL_c` with its quartiles`quantile()` and subset the `consumers` vector correspondingly, and then stored the sampled consumer species' indices in `slct_rand`. Finally, I calculated extinctions and shannon indices similarly as in chapter 2.4 and 2.5 and stored in the respective arrays `extinctions.rand_matTL` and `shannon.rand_matT`.


``` r
shannon.rand_matTL      <- array(NA, dim = c( length(con), length(sample_size) )) 
extinctions.rand_matTL  <- array(NA, dim = c( length(con), length(sample_size) ))

for (j  in 1:length(con)) {
  bas <- n_tot - sum(troph.lvl[j,] > 1)       # the number of basal species
  consumers <- c((bas+1):n_tot)               # a vector of all consumers
  TL_c <- troph.lvl[j,consumers]              # trophic level of each consumer 
  
  for (k in 1:length(sample_size)) {
    
    slct_rand <-  sample(consumers[TL_c <= quantile(TL_c, 0.25)],
                          size = sample_size[k]/4)
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.25) &
                                    TL_c <= quantile(TL_c, 0.5)],
                          size = sample_size[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.5) &
                                    TL_c <= quantile(TL_c, 0.75)],
                          size = sample_size[k]/4))
    
    slct_rand <- c(slct_rand, 
                   sample(consumers[TL_c >= quantile(TL_c, 0.75)],
                          size = sample_size[k]/4))
    
    extinctions.rand_matTL[j,k] <- sum(extinctions[j, slct_rand]) 
    shannon.rand_matTL[j,k] <- diversity((abundances[j, slct_rand]))
  }
}
```

### 2.7 Statistical analysis

## Results

sa

## Discussion

##code chunk \`\`\`{r}

\`\`\`\`r

\`\`\`
