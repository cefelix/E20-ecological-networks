# Title:How does a selection of species represent parameters of the entire food web based on selection criteria?

## Abstract

## Introduction

Many ecological parameters are not measured directly, but assessed through indicator species.

Food webs blabla ...

## 2 Methods

I did all data generation and analysis using `R 4.2.2` (R Core Team, 2022). Food web dynamics were simulated using the `ATNr` package by Gauzens et al. (2022). Other packages used are `vegan` (Oksanen et al., 2022), `ggplot2` (Wickham, 2016), and `tidyr` (Wickham et al., 2023). The complete code for this project can be found at [github.com/cefelix/e20-ecological-networks](https://github.com/cefelix/E20-ecological-networks).

### 2.1 Simulation of food web dynamics

I simulated food web dynamics for 2000 food webs, each food web consisting of `n_tot = 96` trophic species. Connectance `con` varied from 0.05 to 0.35, which corresponds roughly to the span of connectance in empirical food webs as reported by Dunne et al. (2002). Also, I set up the variables `times`, `biom`, and `ext_tresh` which are explained later in this chapter.

``` r
n_tot = 96                                          #total number of species
con <- runif(2000, min=0.05, max=0.35)  %>% sort()  #2000 connectance values

times <- seq(1, 1e12, by = 1e9)   #time for integration of biomass dynamics
biom <- runif(n_tot, 1, 4)        #initial biomasses
ext_thresh <- 0.1**6              #threshold below which considered extinct
```

Next, I generated an adjacency matrix `fw` generated using `ATNr::create_niche_model()`. I assessed the number of basal species `n_bas` by calculating `colSums(fw)`. To minimize effects of variation in basal species number and thus disentangle effects of variation in connectance from effects of variation in basal species number, adjacency matrices were restricted to have between 4 and 12 basal species. This arbitrary range helped avoiding extreme numbers of basal species (the range without restriction lay between 1 and 20+ basal species, with food webs of low connectance typically having very few basal species, while highly connected ones having plenty). Then I set up a body mass vector `BM` with each element representing the body mass of one species. `BM` spanned 12 orders of magnitude, which is considered realistic in soil food webs (Potapov, 2021). Finally, multiple arrays to store output variables were set up prior using `array(dim = c(length(con), n_tot))`. All following steps in chapters 2.1 and 2.2 are to be carried out within the loop through `con` (shown in the code section below).

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

### 2.2 Calculation of output variables

I extracted the final biomasses `bioms` for each species using `sol[nrow(sol), -1]` and stored them in `biomass_array[i, ]` (the first dimension of `biomass_array` corresponds to connectance, while the second corresponds to the species). I applied the same principle to store all the other output variables.

``` r
bioms = sol[nrow(sol),-1]       #final biomasses for each species
biomass_array[i,] = bioms       #storing in array at i-th element of connectance
```

Abundances were calculated dividing `bioms` by individual body masses: `abuns = bioms / BM`. Extinction of each single species was assessed by comparing `bioms` to the previously set extinction threshold 1e-6 using `exts = ifelse(bioms <= model$ext, yes = 1, no = 0)`. Trophic level for each species was calculated with `TroLev(fw)` according to the general formula *INSERT*. The number of initial trophic links (number of prey species `prey.on_start` and consumer species `consumed.by_start`) for each species was calculated by `colSums(fw)` and `rowsums(fw)`, respectively.

``` r
abuns = bioms / BM                                  #final abundances
exts = ifelse(bioms <= model$ext, yes = 1, no = 0)  #extinctions
trolev =  TroLev(fw)                                #trophic levels
prey.on_start = colSums(fw)                         #trophic links: predator
consumed.by_start = rowSums(fw)                     #trophic links: consumed 
```

I then cleared the initial adjacency matrix `fw` from the interactions that involved at least one extinct species. Therefore, I set rows and columns of `fw` to zero in case the corresponding species had become extinct, and stored the result in `fw_end`. The number of final trophic links were calculated as `prey.on_end = colSums(fw_end)` and `consumed.by_end = rowSums(fw_end)`.

``` r
fw_end <- fw
fw_end[exts == 1,] <- 0       #if extinct species: row is set to zero
fw_end[,exts == 1] <- 0       #if extinct species: column is set to zero

prey.on_end = colSums(fw_end)               #trophic links: predator          
consumed.by_end = rowSums(fw_end)           #trophic links: consumed
```

## Results

## Discussion

##code chunk \`\`\`{r}

\`\`\`\`r

\`\`\`
