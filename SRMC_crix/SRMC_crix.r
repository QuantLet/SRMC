#' @usage  Reform the data and save it as CRIX data frame: date, price used for the empirical study
#' @input  crix.json is already saved in the working directory (availabe at 'http://crix.hu-berlin.de')
#' @output SRMC_crix.csv, the log returns of CRIX data.frame with date and price information
#' 
#' 
# clear variables and install necessary packages
# set working directory by setwd('C:/Users/...')
# clear variables and install package
rm(list = ls(all = TRUE))
library("rjson")


crix      = fromJSON(file = "SRMC_crix.json")
size.crix = length(crix)
n         = 1:size.crix

price.crix = function(x, crix) 
    price  = crix[[x]]$price
date.crix  = function(x, crix) 
  date     = crix[[x]]$date

crix.price = sapply(n, price.crix, crix)
crix.date  = sapply(n, date.crix, crix)
# Get the sequence of price and date of the crix

price.log.return = diff(log(crix.price))
# Get the log return of he crix.price

CRIX = data.frame(date = crix.date[2:size.crix], price = price.log.return)

write.csv(CRIX, file = "SRMC_crix", row.names = FALSE) 
# Save the log returns of CRIX as SRMC_crix.csv in the current working directory
