## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE-----------------------------------------------------------
# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)# library(CVXR)
library(pracma)
library(quadprog)
library(parallel)
library(stats)
library(data.table)
library(CVXR)
library(DiSCo)
# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)

## -----------------------------------------------------------------------------
data("dube")
head(dube)

## -----------------------------------------------------------------------------
id_col.target <- 2
t0 <- 2003

## ---- cache=TRUE--------------------------------------------------------------
# if the below gives issues it's the knit cache...
results <- DiSCo(dube, id_col.target, t0, M = 1000, G = 1000, num.cores = 5, permutation = FALSE,
                 CI = TRUE, boots = 1000, cl = 0.95, CI_placebo=TRUE, graph = TRUE, qmethod=NULL)

