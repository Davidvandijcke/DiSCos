# Code for replicating the images in section 5 of the paper "Distributional Synthetic Controls"
# by Florian Gunsilius, Siyun He, and David Van Dijcke, February 2023.

# This code uses publicly available data provided as a supplement to the published article
# "Minimum Wages and the Distribution of Family Incomes" (Dube, 2019, AEJ: Applied).
# The replication package to that article along with the data can be accessed at
# https://www.dropbox.com/s/hrmz2fu0lxet97c/Minimum%20Wage%20Poverty%20Replication.zip?dl=0
# We converted the .dta file to compressed csv to reduce the filesize

# This code comes with no guarantees. The compressed package also contains an unused function
# DSC_CI, which is not used in this replication, but can be used to compute the confidence
# intervals of the corresponding quantile functions.


#####
# Load required packages
packages_load <- c("haven", "base", "data.table", "latex2exp", "CVXR",  # used to compute the weights using the alternative using mixtures of CDF
                   "here", "dplyr", "pracma", "quadprog", "R.utils")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packages_load, character.only = TRUE)


# Set the working directory
setwd(file.path(here::here()))



#####
# Values to be chosen by the researcher
M <- 1000 # draws for samples to compute the respective integrals.
fips.target <- 2 # target state is AK


#####
# Required functions
# loading the function for computing the optimal weights
source('R/DSC_weights_reg.R')
# loading the function for computing the barycenter and the donor distributions
source('R/DSC_bc.R')
# loading the function for performing the permutation test
source('R/DSC_per.R')
## function to compute quantile function
myquant <- function(X,q){
  # sort if unsorted
  if (is.unsorted(X)) X <- sort(X)
  # compute empirical CDF
  X.cdf <- 1:length(X) / length(X)
  # obtain the corresponding empirical quantile
  return(X[which(X.cdf >= q)[1]])
}

#####
# loading the data-set on minimum wage from Dube (2019) to obtain states that did not have a
# change of the minimum wage between 1998-2014
df <- read_dta(file.path("data", "mw_annual_lagsleads_1974_2014.dta"))

years.aff <- list()
for (ii in unique(df$state_fips)){
  hh <- df$mw[df$state_fips==ii & df$year>=1996]
  hh.diff <- diff(hh)
  # getting the years for which the difference is 0
  years.help <- df$year[df$state_fips==ii & df$year>=1996]
  years.aff[[ii]] <- years.help[which(hh.diff==0)]
}

# AK is our treated unit. We want states that did not have a change between 1998 and 2003/2004
states.aff <- list()
for (ii in 1:length(years.aff)){
  dummyvec <- c(1998,1999,2000,2001,2002,2003,2004)
  test <- intersect(dummyvec,years.aff[[ii]])
  states.aff[[ii]] <- FALSE
  if (length(dummyvec) == length(test)){
    if (all(dummyvec == test) == TRUE) {
      states.aff[[ii]] <- TRUE
    }
  }
}
control.states <- which(states.aff == TRUE)

# dropping the data-set for obtaining the control states and loading the actual data
rm(df, states.aff,years.aff, dummyvec,hh,hh.diff,ii,test,years.help)

reread_dta = FALSE # redownload the replication data from the Dropbox package
if (reread_dta) {

  # download the replication data
  url <- "https://www.dropbox.com/s/hrmz2fu0lxet97c/Minimum%20Wage%20Poverty%20Replication.zip?dl=1"
  fn <- file.path("data", "march_regready_1996.zip")
  download.file(url)
  unzip(fn, file.path("Minimum Wage Poverty Replication", "Data", "march_regready_1996.dta"), exdir = "data", junkpaths = TRUE)
  file.remove(fn)

  # load the data from the stata file
  fn <- file.path("data", "march_regready_1996.dta")
  df <- read_dta(fn) %>% setDT()
  df <- df[, c("age", "educ", "contpov", "adj0contpov", "adj1contpov", "adj2contpov", "state_fips", "year", "hhseq", "demgroup1")]
  fwrite(df, file.path("data", "march_regready_1996.csv.gz"))

  file.remove(fn) # remove file to avoid github issues

} else{
  df <- fread(file.path("data", "march_regready_1996.csv.gz"))
}

# only considering individuals under the age of 65
df <- df[demgroup1 == 1]



# collapsing the data
df[,.(age=mean(age),
      educ=mean(educ),
      contpov=mean(contpov),
      adj0contpov=mean(adj0contpov),
      adj1contpov=mean(adj1contpov),
      adj2contpov=mean(adj2contpov)),
   by=c('state_fips', 'year', 'hhseq')]


results.over.years <- list()

#####
# Solving for the optimal weights in the DSC method and the alternative method using mixtures of
# distributions by year

# looping over the years from 1998 - 2014
for(yy in 1:7){
  # obtaining the target state
  target <- list(3)
  target[[3]] <- df[state_fips==fips.target & year==(1997+yy)]$adj0contpov # for the first outcome,
  # adj0contpov
  v <- as.vector(c(rep(TRUE, floor(0.5*length(target[[3]]))),
                   rep(FALSE,ceiling(0.5*length(target[[3]])))))
  set.seed(1860)
  ind <- sample(v)
  target[[1]] <- target[[3]][ind]
  target[[2]] <- target[[3]][!ind]



  # obtaining all control states
  controlstates <- df[state_fips%in%control.states & year==(1997+yy)]
  controls <- list()
  helpcont <- unique(controlstates$state_fips)
  for (ii in 1:length(helpcont)){
    controls[[ii]] <- list(3)
    controls[[ii]][[3]] <- controlstates[state_fips==helpcont[ii] & year==(1997+yy)]$adj0contpov
    # for the first outcome,
    # use adj0contpov
    v <- as.vector(c(rep(TRUE, floor(0.5*length(controls[[ii]][[3]]))),
                     rep(FALSE,ceiling(0.5*length(controls[[ii]][[3]])))))
    set.seed(1860)
    ind <- sample(v)
    controls[[ii]][[1]] <- controls[[ii]][[3]][ind]
    controls[[ii]][[2]] <- controls[[ii]][[3]][!ind]
  }

  # generating the target
  target.s <- mapply(myquant, seq(from=0, to=1, length.out=1001), MoreArgs =
                       list(X=target[[3]]))

  # generating the controls
  controls1 <- list(length(controls))
  for (ii in 1:length(controls)){
    controls1[[ii]] <- controls[[ii]][[3]]
  }


  # obtaining the optimal weights for the DSC method
  DSC_res_weights <- DSC_weights_reg(controls1,as.vector(target[[3]]), 1000)

  DSC_res2 <- DSC_bc(controls1,DSC_res_weights,seq(from=0,to=1,length.out=1001))



  ###### The mixture of distributions approach
  # we again only focus on the first half of the data
  # defining the grid on which we define the cumulative distribution functions
  # obtaining the minimal and maximal values among all supports
  # creating a list of controls with only the full data
  grid.min <- min(c(min(target[[3]]),unlist(lapply(controls1,min))))
  grid.max <- max(c(max(target[[3]]),unlist(lapply(controls1,max))))
  # if grid.max-grid.min<=1 we round to the next 0.1.
  grid.min <- floor(grid.min*10)/10
  grid.max <- ceiling(grid.max*10)/10

  # Estimating the empirical CDFs
  CDF.control <- lapply(controls1,ecdf)
  CDF.target <- ecdf(target[[3]])


  # sampling uniformly on the grid
  grid.rand <- runif(M,grid.min-0.25,grid.max+0.25)
  #grid.rand <- runif(M,0,grid.max+0.25)

  # ordered grid
  grid.ord <- grid.rand[order(grid.rand)]

  # getting the CDF from the quantile function
  DSC_res2.cdfF <- ecdf(DSC_res2[[2]])
  DSC_res2.cdf <- DSC_res2.cdfF(grid.ord)

  # Evaluating the CDF on the random grid
  CDF.matrix <- matrix(0,nrow=length(grid.rand), ncol = (length(controls1)+1))
  CDF.matrix[,1] <- CDF.target(grid.rand)
  for (ii in 1:length(controls1)){
    CDF.matrix[,(ii+1)] <- CDF.control[[ii]](grid.rand)
  }

  # Solving the convex problem with CVXR
  # the variable we are after
  theweights <- Variable(length(controls1))
  # the objective function
  objective <- cvxr_norm((CDF.matrix[,2:ncol(CDF.matrix)] %*% theweights - CDF.matrix[,1]))

  # the constraints for the unit simplex
  constraints <- list(theweights>=0, sum_entries(theweights) == 1)
  # the optimization problem
  problem <- Problem(Minimize(objective),constraints)
  # solving the optimization problem
  results <- solve(problem, solver = "SCS")

  # returning the optimal weights and the value function which provides the
  # squared Wasserstein distance between the target and the corresponding barycenter
  theweights.opt <- results$getValue(theweights)
  thedistance.opt <- results$value*1/M*(grid.max - grid.min)


  themean <- CDF.matrix[,2:ncol(CDF.matrix)]%*%theweights.opt


  themean.order <- themean[order(grid.rand, decreasing=FALSE)]
  target.order <- CDF.matrix[order(grid.rand, decreasing=FALSE),1]

  results.over.years[[yy]] <- list()
  results.over.years[[yy]][[1]] <- list(DSC_res_weights, DSC_res2, DSC_res2.cdf) # DSC estimator
  results.over.years[[yy]][[2]] <- list(theweights.opt, themean, themean.order) # the mixture
  results.over.years[[yy]][[3]] <- list(target.s, target.order, grid.ord, as.vector(target[[3]]))
  results.over.years[[yy]][[4]] <- list(controls1, CDF.matrix)
}

#####
#obtaining the weights as a uniform mean over all time periods
Weights_DSC_avg <- results.over.years[[1]][[1]][[1]]
Weights_mixture_avg <- results.over.years[[1]][[2]][[1]]
for (yy in 2:7){
  Weights_DSC_avg <- Weights_DSC_avg + results.over.years[[yy]][[1]][[1]]
  Weights_mixture_avg <- Weights_mixture_avg + results.over.years[[yy]][[2]][[1]]
}
Weights_DSC_avg <- (1/length(1:7)) * Weights_DSC_avg
Weights_mixture_avg <- (1/length(1:7)) * Weights_mixture_avg

year.to.plot <- 6
# generating the counterfactuals
DSC_res2 <- DSC_bc(results.over.years[[year.to.plot]][[4]][[1]],
                   Weights_DSC_avg,seq(from=0,to=1,length.out=1001))

# getting the CDF from the quantile function
DSC_res2.cdfF <- ecdf(DSC_res2[[2]])
DSC_res2.cdf <- DSC_res2.cdfF(results.over.years[[year.to.plot]][[3]][[3]])

#####
# Plotting the results
pdf("Univ_emp_avg_new_appl.pdf")
plot(results.over.years[[year.to.plot]][[3]][[3]], results.over.years[[year.to.plot]][[2]][[3]],
     type='l', lwd=4,col='#0066FF', xlab='x',ylab='F(x)',cex.lab=1.4, cex.axis=1.4, ylim = c(0,1))
lines(results.over.years[[year.to.plot]][[3]][[3]],DSC_res2.cdf,lwd=4, col='#FF0066', lty=2)
lines(results.over.years[[year.to.plot]][[3]][[3]],results.over.years[[year.to.plot]][[3]][[2]],
      lwd=3, lty=3)
dev.off()

#####
# implementing the permutation test
# setting up the controls over time
controls.per <- list()
target.per <- list()
for (ii in 1:length(results.over.years)){
  controls.per[[ii]] <- results.over.years[[ii]][[4]][[1]]
  target.per[[ii]] <- results.over.years[[yy]][[3]][[4]] # DVD: I think this might be wrong
}


permutation.test.results <- DSC_per(controls.per, target.per, 5,
                                    evgrid=seq(from=0, to=1, length.out=1001), graph=TRUE, y_name='y', x_name='x')

#recording the plot
distt <- permutation.test.results[[1]]
distp <- permutation.test.results[[2]]

#####
# Plot for the permutation test
pdf("Univ_permutation_test_new_appl.pdf")
plot(distt, xaxt ="n", xlab='',ylab='', type='l', lwd=2, ylim = c(0,5))
for (i in 1:length(distp)){
  lines(1:length(distp[[i]]), distp[[i]], col='grey', lwd=1)
}
abline(v=5, lty = 2)
legend("topleft",legend = c("Target", "Control"),
       col=c("black", "grey"),
       lty= c(1,1), lwd = c(2,2), cex = 1.5)
title(ylab="Squared Wasserstein distance", line=2.5, cex.lab=1.5)
title(xlab="Time periods", line=3, cex.lab=1.5)
axis(1, at=1:7, labels=seq(1998,2004,1))
dev.off()








