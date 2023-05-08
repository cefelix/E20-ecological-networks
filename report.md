# Title:How does a selection of species represent parameters of the entire food web based on selection criteria?

## Abstract

## Introduction

Many ecological parameters are not measured directly, but assessed through indicator species.

Food webs blabla ...

## 2 Methods

I did all data generation and analysis using `R 4.2.2` (R Core Team, 2022). Food web dynamics were simulated using the `ATNr` package by Gauzens et al. (2022). Other packages used are `vegan` (Oksanen et al., 2022), `ggplot2` (Wickham, 2016), and `tidyr` (Wickham et al., 2023). The complete code for this project can be found under [github.com/cefelix/e20-ecological-networks](https://github.com/cefelix/E20-ecological-networks).

### 2.1 Simulation of food web dynamics

I simulated food web dynamics for 2000 food webs, each food web consisting of `n_tot = 96` trophic species. Connectance `con` varied from 0.05 to 0.35, which corresponds roughly to the span of connectance in empirical food webs as reported by Dunne et al. (2002). Also, I set up the variables `times`, `biom`, and `ext_tresh` which are explained later in this chapter.

``` r
n_tot = 96                                          #total number of species
con <- runif(2000, min=0.05, max=0.35)  %>% sort()  #2000 connectance values

times <- seq(1, 1e12, by = 1e9)   #time for integration of biomass dynamics
biom <- runif(n_tot, 1, 4)        #initial biomasses
ext_thresh <- 0.1**6              #threshold below which considered extinct
```

Next, I generated an adjacency matrix `fw` generated using `ATNr::create_niche_model()`. I assessed the number of basal species `n_bas` by calculating `colSums(fw)`. To disentangle effects of variation in connectance from effects of variation in basal species number, adjacency matrices were restricted to have between 4 and 12 basal species. Then I set up a body mass vector `BM` with each element representing the bodymass of one species. `BM` spanned 12 orders of magnitude, which is considered realistic in soil food webs (Potapov 2021). All following steps in chapters 2.1 and 2.2 are to be carried out within the here shown loop through `con` as well.

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

Next, an allometric trophic model `model` was created using `create_model_Unscaled()`. Parameters for the ordinary differential equation (ODE) system describing biomass dynamics for each species were initialised by `initialise_model_Unscaled()`. The extinction threshold `model$ext` (i.e. the biomass below which a species is considered extinct) was set to 1e-6 manually. ODE's were then solved numerically over a pre-defined time course `times`, starting from the initial biomasses `bioms` by using `lsoda_wrapper()` and stored in `sol`.

``` r
model <- create_model_Unscaled(n_tot, n_bas, BM, fw) %>%
  initialise_default_Unscaled()
model$ext <- ext_thresh

sol = lsoda_wrapper(times, biom, model) 
```

### 2.2 Calculation of output variables

## Results

## Discussion

##code chunk \`\`\`{r}

´´´
