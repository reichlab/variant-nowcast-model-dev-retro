#' all the functions required for fitting models

require(tidyverse)
require(vcdExtra)
require(rstan)
require(splines)
#' Create a multinomial stan model for a given target date and dataset
#'
#' @param data is a dataframe containing the columns observations, location and date.
#' @param stan_file is the path to which stan file to use
#' @param num_seq is the number of seq required in the last 60 days to be modeled, a numeric integer
#' @param target_date target date is the last day of the data to use, needs to be a date object
#' @param num_days is the number of days before the target date to use, a numeric integer
#' @param target_loc allows you to choose what locations to use, a string vector, default is NULL
#' @param iterations,warmup control the number of iterations and warmup iterations stan uses
#' @param B the number of degrees of freedom if using a spline, default is 3
#' @return return a list containing the stan object, the number of locations, the number of clades, and the locations
stan_maker <- function(data,
                       stan_file,
                       num_seq = 1,
                       target_date = Sys.Date(),
                       num_days = 150,
                       target_loc = NULL,
                       interations = 3000,
                       warmup = 1000,
                       B = 3
                       ){
  # finding which locations to model.
  if( is.null(target_loc)){
    target_lo <- loc_finder(data = data, num_seq = num_seq, target_date = target_date)
  } else{
    target_lo = target_loc
  }
  # filtering to only the locations and dates that we want
  data_case <- filter(data, location %in% target_lo, date >= as.Date(target_date) - num_days, date <= as.Date(target_date))
  # need numeric levels for the mlr
  data_case$mlr <- rep(0, length(data_case$clade))
  # giving each clade a numeric level
  j = 1
  for( k in unique(data$clade)){
    data_case$mlr <- ifelse(data_case$clade == k,j,data_case$mlr )
    j = j + 1
  }
  # need numeric levels for the locations
  data_case$ll <- rep(0, length(data_case$clade))
  j = 1
  for( k in target_lo){
    data_case$ll <- ifelse(data_case$location == k,j,data_case$ll )
    j = j + 1
  }
  # days from start of dataset
  data_case$days <- as.numeric(as_date(data_case$date)) - as.numeric(as_date(as.Date(target_date) - num_days))
  # the number of locations modeled
  L = length(unique(data_case$ll))
  # the number of clades modeled
  K = length(unique(data_case$clade))
  # putting the clades in the right order.
  clades <- unique(data$clade)
  temp <- clades[1]
  clades <- clades[-1]
  clades[K] <- temp
  # creating data for the stan files, checking if we are using a spline model
  if(substr(stan_file, nchar(stan_file)-7, nchar(stan_file)-7) == "S"){
    x <- t(bs(data_case$days, df = B))
    attributes <- attributes(x)
    mlr_data <- list(
      weights = data_case$observation,
      L = L, # number of locations
      ll = data_case$ll, # where each case was located
      y = data_case$mlr, # the clade of each case
      x = x,# the spline over the days the cases happened (first day of dataset treated as 0)
      N = length(data_case$location), # how many cases we had
      K = K, # number of different clades
      B = B # the number of degrees of freedom of the spline
    )
    # fitting the model
    mlr_fit <- stan(
      file = stan_file,
      data = mlr_data,
      chains = 1,
      warmup = warmup,
      iter = interations,
      refresh = 500
    )
    return(list(mlr_fit = mlr_fit, L = L, K = K, target_lo = target_lo, clades = clades, B = B, knots = attributes$knots ))
  } else{
    mlr_data <- list(
      weights = data_case$observation,
      L = L, # number of locations
      ll = data_case$ll, # where each case was located
      y = data_case$mlr, # the clade of each case
      x = data_case$days, # the days the cases happened (first day of dataset treated as 0)
      N = length(data_case$location), # how many cases we had
      K = K # number of different clades
    )
    # fitting the model
    mlr_fit <- stan(
      file = stan_file,
      data = mlr_data,
      chains = 1,
      warmup = warmup,
      iter = interations,
      refresh = 500
    )
    return(list(mlr_fit = mlr_fit, L = L, K = K, target_lo = target_lo, clades = clades ))
  }
}
#' helper function to find which locations match modeling criteria
#'
#' @param data the data file containing the locations
#' @param num_seq the number of sequences required to be modeled
#' @param target_date the date to filter from
#' @returns a vector of strings containing the locations to model
loc_finder <- function(data, num_seq = 1, target_date = Sys.Date()){
  target_lo <- c()
  i = 1
  # finding locations with over num_seq in last 60 days
  for(place in unique(data$location)){
    Virus_counts_state <- filter(data, date >= as.Date(target_date) - 60, location == place)
    if(sum(Virus_counts_state$observation) >= num_seq ){
      target_lo[i] <- place
      i = i + 1
    }
  }
  return(target_lo)
}
#' A helper function to make the data contain only the clades we want for the week
#'
#' @param data the data to transform
#' @param the clades we want to model

trim_clades <- function(data,clades ){
  # seeing if need to create an other clade
  if(any(clades == "other")){
    minor_clades <- c()
    j = 1
    # finding all the clades we don't want to model
    for(clade in unique(data$clade) ){
      if(!(any(clades == clade ))){
        minor_clades[j] <- clade
        j = j + 1
      }
    }
    # merging the clades we don't want to model into the "other" clade
    data$clade <- fct_collapse(data$clade, other = minor_clades)
  } else{
    data <- filter(data, clade %in% clades)
  }
  return(data)
}

#' Create a dirichlet-multinomial stan model for a given target date and dataset
#'
#' @param data is a dataframe containing the columns observations, location and date.
#' @param stan_file is the path to which stan file to use
#' @param num_seq is the number of seq required in the last 60 days to be modeled, a numeric integer
#' @param target_date target date is the last day of the data to use, needs to be a date object
#' @param num_days is the number of days before the target date to use, a numeric integer
#' @param target_loc allows you to choose what locations to use, a string vector, default is NULL
#' @param iterations,warmup control the number of iterations and warmup iterations stan uses
#' @param B the number of degrees of freedom if using a spline, default is 3
#' @return return a list containing the stan object, the number of locations, the number of clades, and the locations
stan_maker_dirichlet <- function(data,stan_file,
                                 num_seq = 1,
                                 target_date = Sys.Date(),
                                 num_days = 150,
                                 target_loc = NULL,
                                 interations = 3000,
                                 warmup = 1000,
                                 B = 3){
  # finding which locations to model.
  if( is.null(target_loc)){
    target_lo <- loc_finder(data = data, num_seq = num_seq, target_date = target_date)
  } else{
    target_lo = target_loc
  }
  # filtering to only the locations and dates that we want
  data_case <- filter(data, location %in% target_lo, date >= as.Date(target_date) - num_days, date <= as.Date(target_date))
  # days from start of dataset
  data_case$days <- as.numeric(as_date(data_case$date)) - as.numeric(as_date(as.Date(target_date) - num_days))
  # the number of locations modeled
  L = length(unique(data_case$location))
  # the number of clades modeled
  clades <- unique(data$clade)
  K = length(clades)
  # creating the array for input
  data_set <- array(dim = c(num_days, L, K))
  # filling in the input array
  for(l in 1:L){
    temp_data <- filter(data_case, location == target_lo[l])
    for(k in 1:K){
      for(n in 1:num_days){
        data_set[n,l,k] <- sum(filter(temp_data, clade == clades[k], days == n )$observation)
      }
    }
  }
  # putting the clades in the right order.
  temp <- clades[1]
  clades <- clades[-1]
  clades[K] <- temp
  # creating data for the stan file
  if(substr(stan_file, nchar(stan_file)-7, nchar(stan_file)-7) == "S"){
    x <- t(bs(1:num_days, df = B))
    attributes <- attributes(x)
  mlr_data <- list(
    N = num_days, # the number of days
    L = L, # the number of locations
    K = K, # the number of clades
    y = data_set, # the counts for each day, location, clade trio
    x = x, # the spline over the days
    B = B # the number of degrees of freedom
  )
  # fitting the model
  mlr_fit <- stan(
    file = stan_file,
    data = mlr_data,
    chains = 1,
    warmup = warmup,
    iter = interations,
    refresh = 500
  )
  return(list(mlr_fit = mlr_fit, L = L, K = K, target_lo = target_lo, clades = clades,B = B, knots = attributes$knots))
  } else{
    mlr_data <- list(
      N = num_days, # the number of days
      L = L, # the number of locations
      K = K, # the number of clades
      y = data_set # the counts for each day, location, clade trio
    )
    # fitting the model
    mlr_fit <- stan(
      file = stan_file,
      data = mlr_data,
      chains = 1,
      warmup = warmup,
      iter = interations,
      refresh = 500
    )
  }
}
