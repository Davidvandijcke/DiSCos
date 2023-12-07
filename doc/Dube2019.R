## ----setup, include = FALSE---------------------------------------------------
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

## ----include=FALSE------------------------------------------------------------

# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)# library(CVXR)
library(pracma)
library(quadprog)
library(parallel)
library(stats)
library(data.table)
library(CVXR)
library(DiSCo)
library(ggplot2)
# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)

## -----------------------------------------------------------------------------
data("dube")
head(dube)

## -----------------------------------------------------------------------------
id_col.target <- 2
t0 <- 2003

## -----------------------------------------------------------------------------
# if the below gives issues it's the knit cache...

disco <- DiSCo(dube, id_col.target, t0, M = 500, G = 1000, num.cores = 2, permutation = FALSE,
                 CI = FALSE, boots = 1, cl = 0.95, CI_placebo=TRUE, graph = FALSE, qmethod=NULL)

## -----------------------------------------------------------------------------
discot <- DiSCoTEA(disco,  agg="quantileDiff", graph=TRUE)
summary(discot)


