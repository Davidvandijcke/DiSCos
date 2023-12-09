## ----setup, include = FALSE----------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
    comment = "#",
    echo=TRUE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE,
  collapse = TRUE,
  out.width = '100%',
  dpi = 144
)


## ----include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------

# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)# library(CVXR)
library(pracma)
library(quadprog)
library(parallel)
library(stats)
library(maps)
library(data.table)
library(CVXR)
library(DiSCo)
library(ggplot2)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
data("dube")
head(dube)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
id_col.target <- 2
t0 <- 2003


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1860)
df <- copy(dube)
disco <- DiSCo(df, id_col.target, t0, G = 1000, num.cores = 5, permutation = FALSE, CI = FALSE, boots = 2, graph = TRUE, simplex=TRUE)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# retrieve the weights
weights <- disco$Weights_DiSCo_avg

# retrieve the control unit IDs
controls <- disco$control_ids

# store in a dataframe
weights_df <- data.frame(weights = weights, fips = controls)

# merge with state fips codes (built into the maps package)
weights_df <- merge(weights_df, maps::state.fips, by = "fips")

setorder(weights_df, -weights)

print(weights_df)



## ------------------------------------------------------------------------------------------------------------------------------------------------------------
maps::state.fips


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(disco$perm)


## ----fig.width=5,fig.height=8, fig.align='center'------------------------------------------------------------------------------------------------------------
discot <- DiSCoTEA(disco,  agg="quantileDiff", graph=TRUE)
summary(discot)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
disco <- DiSCo(dube, id_col.target, t0, G = 200, num.cores = 2, permutation = TRUE, CI = TRUE, boots = 1000, graph = FALSE, q_min = 0.7, q_max=0.95, seed=4)


## ----fig.width=5,fig.height=8, fig.align='center'------------------------------------------------------------------------------------------------------------
discot <- DiSCoTEA(disco, agg="quantileDiff", graph=TRUE)
summary(discot)


