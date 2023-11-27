# Load required packages
packages_load <- c("haven", "base", "data.table",
                   "here", "dplyr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packages_load, character.only = TRUE)


# Set the working directory
setwd(file.path(here::here()))


#####
# loading the data-set on minimum wage from Dube (2019) to obtain states that did not have a
# change of the minimum wage between 1998-2014
url <- "https://www.dropbox.com/s/hrmz2fu0lxet97c/Minimum%20Wage%20Poverty%20Replication.zip?dl=1"
fn <- file.path("data", "mw_annual_lagsleads_1974_2014.zip")
download.file(url, fn)
unzip(fn, file.path("Minimum Wage Poverty Replication", "Data","mw_annual_lagsleads_1974_2014.dta"), exdir = "data", junkpaths = TRUE)
file.remove(fn)
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

# download the replication data
url <- "https://www.dropbox.com/s/hrmz2fu0lxet97c/Minimum%20Wage%20Poverty%20Replication.zip?dl=1"
fn <- file.path("data", "march_regready_1996.zip")
download.file(url, fn)
unzip(fn, file.path("Minimum Wage Poverty Replication", "Data", "march_regready_1996.dta"), exdir = "data", junkpaths = TRUE)
file.remove(fn)

# load the data from the stata file
fn <- file.path("data", "march_regready_1996.dta")
df <- read_dta(fn)  %>% setDT()
# only considering individuals under the age of 65
df <- df[demgroup1 == 1]
df <- df[year %between% c(1998, 2004)]

# collapsing the data
df[,.(adj0contpov=mean(adj0contpov)),
   by=c('state_fips', 'year', 'hhseq')]
dube <- df[, c("adj0contpov", "state_fips", "year")]

fips.target <- 2
dube <- dube[state_fips %in% c(fips.target, control.states)]
dube[, id_col := state_fips]

# create t_col starting at 1 and increasing by 1 for each year
dube[, time_col := year]
dube[, y_col := adj0contpov]

dube <- dube[, c("time_col", "id_col", "y_col")]

file.remove(fn) # remove file to avoid github issues

usethis::use_data(dube, overwrite=TRUE)
